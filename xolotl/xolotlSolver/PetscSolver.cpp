// Includes
#include "PetscSolver.h"
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>

using namespace xolotlCore;

/*
 C_t =  -D*C_xx + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.

 D*C_xx  - diffusion of He[1-5] and V[1] and I[1]
 F(C)  -   forcing function; He being created.
 R(C)  -   reaction terms   (clusters combining)
 D(C)  -   dissociation terms (cluster breaking up)

 Sample Options:
 -ts_monitor_draw_solution               -- plot the solution for each concentration as a function of x each in a separate 1d graph
 -draw_fields_by_name 1-He-2-V,1-He  -- only plot the solution for these two concentrations
 -mymonitor                              -- plot the concentrations of He and V as a function of x and cluster size (2d contour plot)
 -da_refine <n=1,2,...>                  -- run on a finer grid
 -ts_max_steps maxsteps                  -- maximum number of time-steps to take
 -ts_final_time time                     -- maximum time to compute to

 Rules for maximum number of He allowed for V in cluster


 */

namespace xolotlSolver {

//! Help message
static char help[] =
		"Solves C_t =  -D*C_xx + F(C) + R(C) + D(C) from Brian Wirth's SciDAC project.\n";

// Allocate the static network
std::shared_ptr<PSIClusterReactionNetwork> PetscSolver::network;

/**    Hard wire the number of cluster sizes for He, V, and I */
#define N 15

/**
 Define all the concentrations (there is one of these unions at each grid point)

 He[He] represents the clusters of pure Helium of size He
 V[V] the Vacencies of size V,
 I[I] represents the clusters of Interstials of size I,  and
 HeV[He][V]  the mixed Helium-Vacancy clusters of size He and V

 The variables He, V, I are always used to index into the concentrations of He, V, and I respectively
 Note that unlike in traditional C code the indices for He[], V[] and I[] run from 1 to N, NOT 0 to N-1
 (the use of the union below "tricks" the C compiler to allow the indices to start at 1.)

 */
typedef struct {
	PetscScalar He[N];
	PetscScalar V[N];
	union {
		PetscScalar I[N];
		PetscScalar HeV[N + 1][N]; /* actual size is N by N, the N+1 is there only to "trick" the compiler to have indices start at 1.*/
	};
} Concentrations;

/**
 * A simple struct for holding concentrations from the network to pass to
 * PETSc. This structure should be used for networks whose reactants can be
 * ordered sequentially in size.
 */
class Concentrations2 {
public:
	// Array of He
	std::vector<PetscScalar> He;
	// Array of V
	std::vector<PetscScalar> V;
	// Array of I
	std::vector<PetscScalar> I;
	// Array of HeV
	std::vector<PetscScalar> HeV;
	// Array of HeI
	std::vector<PetscScalar> HeI;
};

typedef struct {
	PSIClusterReactionNetwork network;
} PSIClusters;

/**
 * This operation converts a PSIClusterReactionNetwork to Concentrations
 * struct. It assumes that the structures are aligned dimensionally.
 * @param The PSIClusterReactionNetwork whose concentrations should be copied
 * to the concentrations
 * @param The concentrations that should be filled
 */
void convertNetworkToConcentrations(
		const std::shared_ptr<PSIClusterReactionNetwork> & network,
		Concentrations2 * concentrations) {

	// Local Declarations
	int size = 0, i = 0, heSize = 0, vSize = 0, iSize = 0;
	std::map<std::string, int> composition;

	// ----- Convert He -----

	// Setup the He concentrations array
	concentrations->He.clear();
	// Convert the He
	auto reactants = network->getAll("He");
	size = reactants->size();
	// Resize the array
	concentrations->He.resize(size, 0);
	// Copy the He concentrations
	for (i = 0; i < size; i++) {
		concentrations->He.at(i) = reactants->at(i)->getConcentration();
	}

	// ----- Convert V -----

	// Setup the V concentrations array
	concentrations->V.clear();
	// Convert the V
	reactants = network->getAll("V");
	size = reactants->size();
	// Resize the V array
	concentrations->V.resize(size, 0);
	// Copy the V concentrations
	for (i = 0; i < size; i++) {
		concentrations->V.at(i) = reactants->at(i)->getConcentration();
	}

	// ----- Convert I -----

	// Setup the I concentrations array
	concentrations->I.clear();
	// Convert the I
	reactants = network->getAll("I");
	size = reactants->size();
	// Resize the I array
	concentrations->I.resize(size, 0);
	// Copy the I concentrations
	for (i = 0; i < size; i++) {
		concentrations->I.at(i) = reactants->at(i)->getConcentration();
	}

	// ----- Convert HeV -----

	// Setup the HeV concentrations array
	concentrations->HeV.clear();
	// Convert the HeV
	reactants = network->getAll("HeV");
	size = reactants->size();
	// Resize the HeV array
	concentrations->HeV.resize(size, 0);
	// Copy the HeV concentrations
	for (i = 0; i < size; i++) {
		concentrations->HeV.at(i) = reactants->at(i)->getConcentration();
	}

	// ----- Convert HeI -----

	// Setup the HeI concentrations array
	concentrations->HeI.clear();
	// Convert the HeI
	reactants = network->getAll("HeI");
	size = reactants->size();
	// Resize the HeI array
	concentrations->HeI.resize(size, 0);
	// Copy the HeI concentrations
	for (i = 0; i < size; i++) {
		concentrations->HeI.at(i) = reactants->at(i)->getConcentration();
	}

	return;
}

/**
 * This operation converts Concentrations to a PSIClusterReactionNetwork. It
 * assumes that the network has previously been filled with clusters and is
 * aligned dimensionally with the concentrations.
 * @param The concentrations that should be converted
 * @param The PSIClusterReactionNetwork that should be filled
 */
void convertConcentrationsToNetwork(Concentrations2 & concentrations,
		const std::shared_ptr<PSIClusterReactionNetwork> & network) {

	// Local Declarations
	int size = 0, i = 0;
	double conc = 0.0;

	// ----- Convert He -----
	size = concentrations.He.size();
	for (int i = 0; i < size; i++) {
		conc = concentrations.He.at(i);
		network->get("He", i)->setConcentration(conc);
	}

	// ----- Convert V -----
	size = concentrations.V.size();
	for (int i = 0; i < size; i++) {
		conc = concentrations.V.at(i);
		network->get("V", i)->setConcentration(conc);
	}

	// ----- Convert I -----
	size = concentrations.I.size();
	for (int i = 0; i < size; i++) {
		conc = concentrations.I.at(i);
		network->get("I", i)->setConcentration(conc);
	}

	// ----- Convert HeV -----
	size = concentrations.HeV.size();
	auto reactants = network->getAll("HeV");
	for (int i = 0; i < size; i++) {
		conc = concentrations.HeV.at(i);
		reactants->at(i)->setConcentration(conc);
	}

	// ----- Convert HeI -----
	size = concentrations.HeI.size();
	reactants = network->getAll("HeI");
	for (int i = 0; i < size; i++) {
		conc = concentrations.HeI.at(i);
		reactants->at(i)->setConcentration(conc);
	}

	return;
}

/**
 Holds problem specific options and data
 */
typedef struct {
	PetscBool noreactions; /* run without the reaction terms */
	PetscBool nodissociations; /* run without the dissociation terms */
	PetscScalar HeDiffusion[6];
	PetscScalar VDiffusion[2];
	PetscScalar IDiffusion[2];
	PetscScalar forcingScale;
	PetscScalar reactionScale;
	PetscScalar dissociationScale;
} AppCtx;

extern PetscErrorCode IFunction(TS, PetscReal, Vec, Vec, Vec, void*);
extern PetscErrorCode MyMonitorSetUp(TS);

TS ts; /* nonlinear solver */
Vec C; /* solution */
PetscErrorCode ierr;
DM da; /* manages the grid data */
AppCtx ctx; /* holds problem specific parameters */
PetscInt He, dof = 3 * N + N * N, *ofill, *dfill;

/* ----- Error Handling Code ----- */

/**
 * This operation checks a Petsc error code and converts it to a bool.
 * @param errorCode The Petsc error code.
 * @return True if everything is OK, false otherwise.
 */
static inline bool checkPetscError(PetscErrorCode errorCode) {
	CHKERRQ(errorCode);
}

/**
 * This operation "returns" in a way that Petsc expects.
 * @return The return code from Petsc.
 */
static inline int petscReturn() {
	PetscFunctionReturn(0);
}

/* ----- Monitoring Code ----- */

typedef struct {
	DM da; /* defines the 2d layout of the He subvector */
	Vec He;
	VecScatter scatter;
	PetscViewer viewer;
} MyMonitorCtx;

#undef __FUNCT__
#define __FUNCT__ "monitorSolve"
/**
 * This is a monitoring operation that displays He and V as a function of space
 * and cluster size for each time step. It is not a member variable of the class
 * because the monitoring code requires a C callback function (via a function
 * pointer).
 */
static PetscErrorCode monitorSolve(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx) {
	MyMonitorCtx *ctx = (MyMonitorCtx*) ictx;
	PetscErrorCode ierr;

	PetscFunctionBeginUser;
	ierr = VecScatterBegin(ctx->scatter, solution, ctx->He, INSERT_VALUES,
			SCATTER_FORWARD);
	checkPetscError(ierr);
	ierr = VecScatterEnd(ctx->scatter, solution, ctx->He, INSERT_VALUES,
			SCATTER_FORWARD);
	checkPetscError(ierr);
	ierr = VecView(ctx->He, ctx->viewer);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "destroyPetscMonitor"
/**
 * This operation frees all data structures associated with the monitor. It is
 * not a member variable of the class because the monitoring code requires a C
 * callback function (via a function pointer).
 * @param ictx The monitor
 * @return A standard PETSc error code
 */
static PetscErrorCode destroyPetscMonitor(void **ictx) {
	MyMonitorCtx **ctx = (MyMonitorCtx**) ictx;
	PetscErrorCode ierr;

	PetscFunctionBeginUser;
	ierr = VecScatterDestroy(&(*ctx)->scatter);
	checkPetscError(ierr);
	ierr = VecDestroy(&(*ctx)->He);
	checkPetscError(ierr);
	ierr = DMDestroy(&(*ctx)->da);
	checkPetscError(ierr);
	ierr = PetscViewerDestroy(&(*ctx)->viewer);
	checkPetscError(ierr);
	ierr = PetscFree(*ctx);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "setupPetscMonitor"
/**
 * This operation sets up a monitor that will display He as a function of space
 * and cluster size for each time step. It was not made a member function so
 * that it would be consistent with the other monitor callbacks.
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
static PetscErrorCode setupPetscMonitor(TS ts) {
	DM da;
	PetscErrorCode ierr;
	PetscInt xi, xs, xm, *idx, M, xj, cnt = 0, dof = 3 * N + N * N;
	const PetscInt *lx;
	Vec C;
	MyMonitorCtx *ctx;
	PetscBool flg;
	IS is;
	char ycoor[32];
	PetscReal valuebounds[4] = { 0, 1.2, 0, 1.2 };

	PetscFunctionBeginUser;
	ierr = PetscOptionsHasName(NULL, "-mymonitor", &flg);
	checkPetscError(ierr);
	if (!flg)
		PetscFunctionReturn(0);

	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);
	ierr = PetscNew(MyMonitorCtx, &ctx);
	checkPetscError(ierr);
	ierr = PetscViewerDrawOpen(PetscObjectComm((PetscObject) da), NULL, "",
			PETSC_DECIDE, PETSC_DECIDE, 600, 400, &ctx->viewer);
	checkPetscError(ierr);

	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &M, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);
	ierr = DMDAGetOwnershipRanges(da, &lx, NULL, NULL);
	checkPetscError(ierr);
	ierr = DMDACreate2d(PetscObjectComm((PetscObject) da), DMDA_BOUNDARY_NONE,
			DMDA_BOUNDARY_NONE, DMDA_STENCIL_STAR, M, N, PETSC_DETERMINE, 1, 2,
			1, lx, NULL, &ctx->da);
	checkPetscError(ierr);
	ierr = DMDASetFieldName(ctx->da, 0, "He");
	checkPetscError(ierr);
	ierr = DMDASetFieldName(ctx->da, 1, "V");
	checkPetscError(ierr);
	ierr = DMDASetCoordinateName(ctx->da, 0, "X coordinate direction");
	checkPetscError(ierr);
	ierr = PetscSNPrintf(ycoor, 32, "%D ... Cluster size ... 1", N);
	checkPetscError(ierr);
	ierr = DMDASetCoordinateName(ctx->da, 1, ycoor);
	checkPetscError(ierr);

	ierr = DMCreateGlobalVector(ctx->da, &ctx->He);
	checkPetscError(ierr);
	ierr = PetscMalloc(2 * N * xm * sizeof(PetscInt), &idx);
	checkPetscError(ierr);
	cnt = 0;
	for (xj = 0; xj < N; xj++) {
		for (xi = xs; xi < xs + xm; xi++) {
			idx[cnt++] = dof * xi + xj;
			idx[cnt++] = dof * xi + xj + N;
		}
	}
	ierr = ISCreateGeneral(PetscObjectComm((PetscObject) ts), 2 * N * xm, idx,
			PETSC_OWN_POINTER, &is);
	checkPetscError(ierr);
	ierr = TSGetSolution(ts, &C);
	checkPetscError(ierr);
	ierr = VecScatterCreate(C, is, ctx->He, NULL, &ctx->scatter);
	checkPetscError(ierr);
	ierr = ISDestroy(&is);
	checkPetscError(ierr);

	/* sets the bounds on the contour plot values so the colors mean the same thing for different timesteps */
	ierr = PetscViewerDrawSetBounds(ctx->viewer, 2, valuebounds);
	checkPetscError(ierr);

	ierr = TSMonitorSet(ts, monitorSolve, ctx, destroyPetscMonitor);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "InitialConditions"
PetscErrorCode PetscSolver::setupInitialConditions(DM da, Vec C) {

	// Local Declarations
	PetscErrorCode ierr;
	PetscInt i, nI, nHe, nV, xs, xm, Mx, cnt = 0;
	Concentrations2 *concentrations;
	PetscReal hx, x;
	char string[16];
	auto reactants = network->getAll();
	int size = reactants->size();

	PetscFunctionBeginUser;
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);
	hx = 1.0 / (PetscReal) (Mx - 1);

	/* Name each of the concentrations */
	for (i = 0; i < size; i++) {
		auto composition = reactants->at(i)->getComposition();
		nHe = composition["He"];
		nV = composition["V"];
		nI = composition["I"];
		ierr = PetscSNPrintf(string, 16, "He-%d,V-%d,I-%d", nHe, nV, nI);
		checkPetscError(ierr);
		ierr = DMDASetFieldName(da, cnt++, string);
		checkPetscError(ierr);
	}

//	for (He = 1; He < N + 1; He++) {
//		ierr = PetscSNPrintf(string, 16, "%d-He", He);
//		checkPetscError(ierr);
//		ierr = DMDASetFieldName(da, cnt++, string);
//		checkPetscError(ierr);
//	}
//	for (V = 1; V < N + 1; V++) {
//		ierr = PetscSNPrintf(string, 16, "%d-V", V);
//		checkPetscError(ierr);
//		ierr = DMDASetFieldName(da, cnt++, string);
//		checkPetscError(ierr);
//	}
//	for (I = 1; I < N + 1; I++) {
//		ierr = PetscSNPrintf(string, 16, "%d-I", I);
//		checkPetscError(ierr);
//		ierr = DMDASetFieldName(da, cnt++, string);
//		checkPetscError(ierr);
//	}
//	for (He = 1; He < N + 1; He++) {
//		for (V = 1; V < N + 1; V++) {
//			ierr = PetscSNPrintf(string, 16, "%d-He-%d-V", He, V);
//			checkPetscError(ierr);
//			ierr = DMDASetFieldName(da, cnt++, string);
//			checkPetscError(ierr);
//		}
//	}

	/*
	 Get pointer to vector data
	 */
	ierr = DMDAVecGetArray(da, C, &concentrations);
	checkPetscError(ierr);

	/*
	 Get local grid boundaries
	 */
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	/*
	 Compute function over the locally owned part of the grid
	 */
	for (i = xs; i < xs + xm; i++) {
		x = i * hx;
		// Set the default vacancy concentrations
		reactants = network->getAll("V");
		size = reactants->size();
		for (int j = 0; j < size; j++) {
			reactants->at(j)->setConcentration(1.0);
		}
		// Set the default interstitial concentrations
		reactants = network->getAll("I");
		size = reactants->size();
		for (int j = 0; j < size; j++) {
			reactants->at(j)->setConcentration(1.0);
		}
		// Create a copy of the network for this grid point by converting to a
		// Concentrations struct
		convertNetworkToConcentrations(network, &concentrations[i]);
//		for (He = 1; He < N + 1; He++)
//			c[i].He[He] = 0.0;
//		for (V = 1; V < N + 1; V++)
//			c[i].V[V] = 1.0;
//		for (I = 1; I < N + 1; I++)
//			c[i].I[I] = 1.0;
//		for (He = 1; He < N + 1; He++) {
//			for (V = 1; V < N + 1; V++)
//				c[i].HeV[He][V] = 0.0;
//		}
	}

	/*
	 Restore vectors
	 */
	concentrations = (Concentrations2*) (((PetscScalar*) concentrations) + 1);
	ierr = DMDAVecRestoreArray(da, C, &concentrations);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}

/* ------------------------------------------------------------------- */
#undef __FUNCT__
#define __FUNCT__ "IFunction"
/*
 IFunction - Evaluates nonlinear function that defines the ODE

 Input Parameters:
 .  ts - the TS context
 .  U - input vector
 .  ptr - optional user-defined context

 Output Parameter:
 .  F - function values
 */
PetscErrorCode IFunction(TS ts, PetscReal ftime, Vec C, Vec Cdot, Vec F,
		void *ptr) {
	AppCtx *ctx = (AppCtx*) ptr;
	DM da;
	PetscErrorCode ierr;
	PetscInt xi, Mx, xs, xm, He, he, V, v, I, i;
	PetscReal hx, sx, x;
	Concentrations2 *concs, *updatedConcs;
	Vec localC;
	std::shared_ptr<PSICluster> newCluster;
	std::shared_ptr<std::vector<std::shared_ptr<Reactant>>>oldReactants, newReactants;
	int size = 0;
	double oldConc = 0.0, oldLeftConc = 0.0, oldRightConc = 0.0, conc = 0.0,
			temperature = 1000.0, flux = 0.0;

	PetscFunctionBeginUser;
	ierr = TSGetDM(ts, &da);
	checkPetscError(ierr);
	ierr = DMGetLocalVector(da, &localC);
	checkPetscError(ierr);
	ierr = DMDAGetInfo(da, PETSC_IGNORE, &Mx, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE, PETSC_IGNORE,
			PETSC_IGNORE);
	checkPetscError(ierr);
	hx = 8.0 / (PetscReal) (Mx - 1);
	sx = 1.0 / (hx * hx);

	/*
	 F  = Cdot +  all the diffusion and reaction terms added below
	 */
	ierr = VecCopy(Cdot, F);
	checkPetscError(ierr);

	/*
	 Scatter ghost points to local vector,using the 2-step process
	 DMGlobalToLocalBegin(),DMGlobalToLocalEnd().
	 By placing code between these two statements, computations can be
	 done while messages are in transition.
	 */
	ierr = DMGlobalToLocalBegin(da, C, INSERT_VALUES, localC);
	checkPetscError(ierr);
	ierr = DMGlobalToLocalEnd(da, C, INSERT_VALUES, localC);
	checkPetscError(ierr);

	/*
	 Get pointers to vector data
	 */
	ierr = DMDAVecGetArray(da, localC, &concs);
	checkPetscError(ierr);
	/* Shift the c pointer to allow accessing with index of 1, instead of 0 */
	ierr = DMDAVecGetArray(da, F, &updatedConcs);
	checkPetscError(ierr);
	updatedConcs = (Concentrations2*) (((PetscScalar*) updatedConcs) - 1);

	/*
	 Get local grid boundaries
	 */
	ierr = DMDAGetCorners(da, &xs, NULL, NULL, &xm, NULL, NULL);
	checkPetscError(ierr);

	// Get the network
	auto network = PetscSolver::getNetwork();

	/*
	 Loop over grid points computing ODE terms for each grid point
	 */
	for (xi = xs; xi < xs + xm; xi++) {
		x = xi * hx;
		// Copy data into the PSIClusterReactionNetwork so that it can
		// compute the new concentrations.
		convertConcentrationsToNetwork(updatedConcs[xi],
				PetscSolver::getNetwork());

		/* ---- Compute diffusion over the locally owned part of the grid -----

		 /* He clusters larger than 5 do not diffuse -- are immobile */
		for (i = 1; i < 6; i++) {
			// Get the concentrations
			oldConc = concs[xi].He.at(i);
			oldLeftConc = concs[xi - 1].He.at(i);
			oldRightConc = concs[xi + 1].He.at(i);
			// Get size*He from the new network
			newCluster = std::dynamic_pointer_cast<PSICluster>(
					network->get("He", size));
			// Only update the concentration if the clusters exist
			if (newCluster) {
				// Use a simple midpoint stencil to compute the concentration
				conc = newCluster->getDiffusionCoefficient(temperature)
						* (-2.0 * oldConc + oldLeftConc + oldRightConc) * sx;
				// Update the concentration of the new cluster
				newCluster->decreaseConcentration(conc);
			}
		}

		/* ----- Compute forcing that produces He of cluster size 1 -----
		 Crude cubic approximation of graph from Tibo's notes
		 */
		concs[xi].He.at(i) -= ctx->forcingScale * PetscMax(0.0,
				0.0006 * x * x * x - 0.0087 * x * x + 0.0300 * x);
		/* Are V or I produced? */

//		for (He = 1; He < PetscMin(N + 1, 6); He++) {
//			f[xi].He[He] -=
//					ctx->HeDiffusion[He]
//							* (-2.0 * c[xi].He[He] + c[xi - 1].He[He]
//									+ c[xi + 1].He[He]) * sx;
//		}
		// Only vacancy clusters of size 1 diffuse.
		size = 1;
		// Get the concentrations.
		oldConc = concs[xi].V.at(i);
		oldLeftConc = concs[xi - 1].V.at(i);
		oldRightConc = concs[xi + 1].V.at(i);
		// Get size*V from the new network
		newCluster = std::dynamic_pointer_cast<PSICluster>(
				network->get("V", size));
		// Only update the concentration if the clusters exist
		if (newCluster) {
			// Use a simple midpoint stencil to compute the concentration
			conc = newCluster->getDiffusionCoefficient(temperature)
					* (-2.0 * oldConc + oldLeftConc + oldRightConc) * sx;
			// Update the concentration of the new cluster
			newCluster->decreaseConcentration(conc);
		}

		// Only interstitial clusters of size 1 diffuse. Get the
		// concentrations.
		oldConc = concs[xi].I.at(i);
		oldLeftConc = concs[xi - 1].I.at(i);
		oldRightConc = concs[xi + 1].I.at(i);
		// Get size*I from the new network
		newCluster = std::dynamic_pointer_cast<PSICluster>(
				network->get("I", size));
		// Only update the concentration if the clusters exist
		if (newCluster) {
			// Use a simple midpoint stencil to compute the concentration
			conc = newCluster->getDiffusionCoefficient(temperature)
					* (-2.0 * oldConc + oldLeftConc + oldRightConc) * sx;
			// Update the concentration of the new cluster
			newCluster->decreaseConcentration(conc);
		}

		/* V and I clusters ONLY of size 1 diffuse */
//		f[xi].V[1] -= ctx->VDiffusion[1]
//				* (-2.0 * c[xi].V[1] + c[xi - 1].V[1] + c[xi + 1].V[1]) * sx;
//		f[xi].I[1] -= ctx->IDiffusion[1]
//				* (-2.0 * c[xi].I[1] + c[xi - 1].I[1] + c[xi + 1].I[1]) * sx;
		/*----- Mixed He - V clusters are immobile  -----*/

//		if (ctx->noreactions)
//			continue;
		// Get all of the old and new reactants
//		oldReactants = clusters[xi].network.getAll();
//		newReactants = updatedClusters[xi].network.getAll();
//		// Update the concentrations based fluxes from reaction term
//		size = oldReactants->size();
//		for (i = 0; i < size; i++) {
//			newReactants->at(i)->decreaseConcentration(
//					oldReactants->at(i)->getTotalFlux(1000));
//		}
		// ----- Convert He -----
		size = concs[xi].He.size();
		for (int i = 0; i < size; i++) {
			newCluster = std::dynamic_pointer_cast<PSICluster>(
					network->get("He", i));
			flux = newCluster->getTotalFlux(temperature);
			newCluster->decreaseConcentration(flux);
		}

		// ----- Convert V -----
		size = concs[xi].V.size();
		for (int i = 0; i < size; i++) {
			newCluster = std::dynamic_pointer_cast<PSICluster>(
					network->get("V", size));
			flux = newCluster->getTotalFlux(temperature);
			newCluster->decreaseConcentration(flux);
		}

		// ----- Convert I -----
		size = concs[xi].I.size();
		for (int i = 0; i < size; i++) {
			newCluster = std::dynamic_pointer_cast<PSICluster>(
					network->get("I", size));
			flux = newCluster->getTotalFlux(temperature);
			newCluster->decreaseConcentration(flux);
		}

		// ----- Convert HeV -----
		size = concs[xi].HeV.size();
		auto reactants = network->getAll("HeV");
		for (int i = 0; i < size; i++) {
			newCluster = std::dynamic_pointer_cast<PSICluster>(
					reactants->at(i));
			flux = newCluster->getTotalFlux(temperature);
			newCluster->decreaseConcentration(flux);
		}

		// ----- Convert HeI -----
		size = concs[xi].HeI.size();
		reactants = network->getAll("HeI");
		for (int i = 0; i < size; i++) {
			newCluster = std::dynamic_pointer_cast<PSICluster>(
					reactants->at(i));
			flux = newCluster->getTotalFlux(temperature);
			newCluster->decreaseConcentration(flux);
			;
		}

		// Convert the concentrations back to the PETSc structure
		convertNetworkToConcentrations(network, &updatedConcs[xi]);

		/* ----------------------------------------------------------------
		 ---- Compute reaction terms that can create a cluster of given size
		 */
		/*   He[He] + He[he] -> He[He+he]  */
//		for (He = 2; He < N + 1; He++) {
//			/* compute all pairs of clusters of smaller size that can combine to create a cluster of size He,
//			 remove the upper half since they are symmetric to the lower half of the pairs. For example
//			 when He = 5 (cluster size 5) the pairs are
//			 1   4
//			 2   2
//			 3   2  these last two are not needed in the sum since they repeat from above
//			 4   1  this is why he < (He/2) + 1            */
//			for (he = 1; he < (He / 2) + 1; he++) {
//				f[xi].He[He] -= ctx->reactionScale * c[xi].He[he]
//						* c[xi].He[He - he];
//
//				/* remove the two clusters that merged to form the larger cluster */
//				f[xi].He[he] += ctx->reactionScale * c[xi].He[he]
//						* c[xi].He[He - he];
//				f[xi].He[He - he] += ctx->reactionScale * c[xi].He[he]
//						* c[xi].He[He - he];
//			}
//		}
//		/*   V[V]  +  V[v] ->  V[V+v]  */
//		for (V = 2; V < N + 1; V++) {
//			for (v = 1; v < (V / 2) + 1; v++) {
//				f[xi].V[V] -= ctx->reactionScale * c[xi].V[v] * c[xi].V[V - v];
//				/* remove the clusters that merged to form the larger cluster */
//				f[xi].V[v] += ctx->reactionScale * c[xi].V[v] * c[xi].V[V - v];
//				f[xi].V[V - v] += ctx->reactionScale * c[xi].V[v]
//						* c[xi].V[V - v];
//			}
//		}
//		/*   I[I] +  I[i] -> I[I+i] */
//		for (I = 2; I < N + 1; I++) {
//			for (i = 1; i < (I / 2) + 1; i++) {
//				f[xi].I[I] -= ctx->reactionScale * c[xi].I[i] * c[xi].I[I - i];
//				/* remove the clusters that merged to form the larger cluster */
//				f[xi].I[i] += ctx->reactionScale * c[xi].I[i] * c[xi].I[I - i];
//				f[xi].I[I - i] += ctx->reactionScale * c[xi].I[i]
//						* c[xi].I[I - i];
//			}
//		}
//		/* He[1] +  V[1]  ->  He[1]-V[1] */
//		f[xi].HeV[1][1] -= 1000 * ctx->reactionScale * c[xi].He[1] * c[xi].V[1];
//		/* remove the He and V  that merged to form the He-V cluster */
//		f[xi].He[1] += 1000 * ctx->reactionScale * c[xi].He[1] * c[xi].V[1];
//		f[xi].V[1] += 1000 * ctx->reactionScale * c[xi].He[1] * c[xi].V[1];
//		/*  He[He]-V[V] + He[he] -> He[He+he]-V[V]  */
//		for (He = 1; He < N; He++) {
//			for (V = 1; V < N + 1; V++) {
//				for (he = 1; he < N - He + 1; he++) {
//					f[xi].HeV[He + he][V] -= ctx->reactionScale
//							* c[xi].HeV[He][V] * c[xi].He[he];
//					/* remove the two clusters that merged to form the larger cluster */
//					f[xi].He[he] += ctx->reactionScale * c[xi].HeV[He][V]
//							* c[xi].He[he];
//					f[xi].HeV[He][V] += ctx->reactionScale * c[xi].HeV[He][V]
//							* c[xi].He[he];
//				}
//			}
//		}
//		/*  He[He]-V[V] + V[v] -> He[He][V+v] */
//		for (He = 1; He < N + 1; He++) {
//			for (V = 1; V < N; V++) {
//				for (v = 1; v < N - V + 1; v++) {
//					f[xi].HeV[He][V + v] -= ctx->reactionScale
//							* c[xi].HeV[He][V] * c[xi].V[v];
//					/* remove the two clusters that merged to form the larger cluster */
//					f[xi].V[v] += ctx->reactionScale * c[xi].HeV[He][V]
//							* c[xi].V[v];
//					f[xi].HeV[He][V] += ctx->reactionScale * c[xi].HeV[He][V]
//							* c[xi].V[v];
//				}
//			}
//		}
//		/*  He[He]-V[V]  + He[he]-V[v] -> He[He+he][V+v]  */
//		/*  Currently the reaction rates for this are zero */
//		for (He = 1; He < N; He++) {
//			for (V = 1; V < N; V++) {
//				for (he = 1; he < N - He + 1; he++) {
//					for (v = 1; v < N - V + 1; v++) {
//						f[xi].HeV[He + he][V + v] -= 0.0 * c[xi].HeV[He][V]
//								* c[xi].HeV[he][v];
//						/* remove the two clusters that merged to form the larger cluster */
//						f[xi].HeV[he][V] += 0.0 * c[xi].HeV[He][V]
//								* c[xi].HeV[he][v];
//						f[xi].HeV[He][V] += 0.0 * c[xi].HeV[He][V]
//								* c[xi].HeV[he][v];
//					}
//				}
//			}
//		}
//		/*  V[V] + I[I]  ->   V[V-I] if V > I else I[I-V] */
//		/*  What should the correct reaction rate should be? */
//		for (V = 1; V < N + 1; V++) {
//			for (I = 1; I < V; I++) {
//				f[xi].V[V - I] -= ctx->reactionScale * c[xi].V[V] * c[xi].I[I];
//				f[xi].V[V] += ctx->reactionScale * c[xi].V[V] * c[xi].I[I];
//				f[xi].I[I] += ctx->reactionScale * c[xi].V[V] * c[xi].I[I];
//			}
//			for (I = V + 1; I < N + 1; I++) {
//				f[xi].I[I - V] -= ctx->reactionScale * c[xi].V[V] * c[xi].I[I];
//				f[xi].V[V] += ctx->reactionScale * c[xi].V[V] * c[xi].I[I];
//				f[xi].I[I] += ctx->reactionScale * c[xi].V[V] * c[xi].I[I];
//			}
//		}
//
//		if (ctx->nodissociations)
//			continue;
//		/* -------------------------------------------------------------------------
//		 ---- Compute dissociation terms that removes an item from a cluster
//		 I assume dissociation means losing only a single item from a cluster
//		 I cannot tell from the notes if clusters can break up into any sub-size.
//		 */
//		/*   He[He] ->  He[He-1] + He[1] */
//		for (He = 2; He < N + 1; He++) {
//			f[xi].He[He - 1] -= ctx->dissociationScale * c[xi].He[He];
//			f[xi].He[1] -= ctx->dissociationScale * c[xi].He[He];
//			f[xi].He[He] += ctx->dissociationScale * c[xi].He[He];
//		}
//		/*   V[V] ->  V[V-1] + V[1] */
//		for (V = 2; V < N + 1; V++) {
//			f[xi].V[V - 1] -= ctx->dissociationScale * c[xi].V[V];
//			f[xi].V[1] -= ctx->dissociationScale * c[xi].V[V];
//			f[xi].V[V] += ctx->dissociationScale * c[xi].V[V];
//		}
//		/*   I[I] ->  I[I-1] + I[1] */
//		for (I = 2; I < N + 1; I++) {
//			f[xi].I[I - 1] -= ctx->dissociationScale * c[xi].I[I];
//			f[xi].I[1] -= ctx->dissociationScale * c[xi].I[I];
//			f[xi].I[I] += ctx->dissociationScale * c[xi].I[I];
//		}
//		/* He[1]-V[1]  ->  He[1] + V[1] */
//		f[xi].He[1] -= 1000 * ctx->reactionScale * c[xi].HeV[1][1];
//		f[xi].V[1] -= 1000 * ctx->reactionScale * c[xi].HeV[1][1];
//		f[xi].HeV[1][1] += 1000 * ctx->reactionScale * c[xi].HeV[1][1];
//		/*   He[He]-V[1] ->  He[He] + V[1]  */
//		for (He = 2; He < N + 1; He++) {
//			f[xi].He[He] -= 1000 * ctx->reactionScale * c[xi].HeV[He][1];
//			f[xi].V[1] -= 1000 * ctx->reactionScale * c[xi].HeV[He][1];
//			f[xi].HeV[He][1] += 1000 * ctx->reactionScale * c[xi].HeV[He][1];
//		}
//		/*   He[1]-V[V] ->  He[1] + V[V]  */
//		for (V = 2; V < N + 1; V++) {
//			f[xi].He[1] -= 1000 * ctx->reactionScale * c[xi].HeV[1][V];
//			f[xi].V[V] -= 1000 * ctx->reactionScale * c[xi].HeV[1][V];
//			f[xi].HeV[1][V] += 1000 * ctx->reactionScale * c[xi].HeV[1][V];
//		}
//		/*   He[He]-V[V] ->  He[He-1]-V[V] + He[1]  */
//		for (He = 2; He < N + 1; He++) {
//			for (V = 2; V < N + 1; V++) {
//				f[xi].He[1] -= 1000 * ctx->reactionScale * c[xi].HeV[He][V];
//				f[xi].HeV[He - 1][V] -= 1000 * ctx->reactionScale
//						* c[xi].HeV[He][V];
//				f[xi].HeV[He][V] += 1000 * ctx->reactionScale
//						* c[xi].HeV[He][V];
//			}
//		}
//		/*   He[He]-V[V] ->  He[He]-V[V-1] + V[1]  */
//		for (He = 2; He < N + 1; He++) {
//			for (V = 2; V < N + 1; V++) {
//				f[xi].V[1] -= 1000 * ctx->reactionScale * c[xi].HeV[He][V];
//				f[xi].HeV[He][V - 1] -= 1000 * ctx->reactionScale
//						* c[xi].HeV[He][V];
//				f[xi].HeV[He][V] += 1000 * ctx->reactionScale
//						* c[xi].HeV[He][V];
//			}
//		}
//		/*   He[He]-V[V] ->  He[He]-V[V+1] + I[1]  */
//		/* Again, what is the reasonable dissociation rate? */
//		for (He = 1; He < N + 1; He++) {
//			for (V = 1; V < N; V++) {
//				f[xi].HeV[He][V + 1] -= 1000 * ctx->reactionScale
//						* c[xi].HeV[He][V];
//				f[xi].I[1] -= 1000 * ctx->reactionScale * c[xi].HeV[He][V];
//				f[xi].HeV[He][V] += 1000 * ctx->reactionScale
//						* c[xi].HeV[He][V];
//			}
//		}
//
	}

	/*
	 Restore vectors
	 */
	ierr = DMDAVecRestoreArray(da, localC, &concs);
	checkPetscError(ierr);
	ierr = DMDAVecRestoreArray(da, F, &updatedConcs);
	checkPetscError(ierr);
	ierr = DMRestoreLocalVector(da, &localC);
	checkPetscError(ierr);
	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "getDiagonalFill"

PetscErrorCode PetscSolver::getDiagonalFill(PetscInt *diagFill,
		int diagFillSize) {

	// Local Declarations
	int i = 0, j = 0, numReactants = network->size(), index = 0,
			connectivityLength = 0, size = numReactants * numReactants;
	std::shared_ptr<std::vector<int> > connectivity;
	std::shared_ptr<Reactant> reactant;

	// Fill the diagonal block if the sizes match up
	if (diagFillSize == size) {
		auto reactants = network->getAll();
		auto testReactants = *reactants;
		// Get the connectivity for each reactant
		for (i = 0; i < numReactants; i++) {
			// Get the reactant and its connectivity
			reactant = reactants->at(i);
			connectivity = reactant->getConnectivity();
			// Add it to the diagonal fill block
			connectivityLength = connectivity->size();
			for (j = 0; j < connectivityLength; j++) {
				index = i * numReactants + j;
				diagFill[index] = connectivity->at(j);
			}
		}
	} else {
		std::string err =
				"PetscSolver Exception: Invalid diagonal block size!\n";
		throw std::string(err);
	}

//
//	if (!ctx->noreactions) {
//
//		for (He = 2; He < N + 1; He++) {
//			/* compute all pairs of clusters of smaller size that can combine to create a cluster of size He,
//			 remove the upper half since they are symmetric to the lower half of the pairs. For example
//			 when He = 5 (cluster size 5) the pairs are
//			 1   4
//			 2   2
//			 3   2  these last two are not needed in the sum since they repeat from above
//			 4   1  this is why he < (He/2) + 1            */
//			for (he = 1; he < (He / 2) + 1; he++) {
//				reactants[0] = he, reactants[1] = He - he, reactants[2] = He;
//				for (j = 0; j < 3; j++) {
//					row = reactants[j], col1 = reactants[0], col2 = reactants[1];
//					int index1 = (row - 1) * dof + col1 - 1, index2 = (row - 1)
//							* dof + col2 - 1;
//					dfill[index1] = 1, dfill[index2] = 1;
//					//std::cout << He << " " << he << " " << row << " " << col1 << " " << col2 << " " << index1 << " " << index2 << std::endl;
//				}
//			}
//		}
//
//		/*   V[V]  +  V[v] ->  V[V+v]  */
//		for (V = 2; V < N + 1; V++) {
//			for (v = 1; v < (V / 2) + 1; v++) {
//				reactants[0] = N + v, reactants[1] = N + V - v, reactants[2] = N
//						+ V;
//				for (j = 0; j < 3; j++) {
//					row = reactants[j], col1 = reactants[0], col2 = reactants[1];
//					dfill[(row - 1) * dof + col1 - 1] = 1, dfill[(row - 1) * dof
//							+ col2 - 1] = 1;
//				}
//			}
//		}
//
//		/*   I[I] +  I[i] -> I[I+i] */
//		for (I = 2; I < N + 1; I++) {
//			for (i = 1; i < (I / 2) + 1; i++) {
//				reactants[0] = 2 * N + i, reactants[1] = 2 * N + I - i, reactants[2] =
//						2 * N + I;
//				for (j = 0; j < 3; j++) {
//					row = reactants[j], col1 = reactants[0], col2 = reactants[1];
//					dfill[(row - 1) * dof + col1 - 1] = 1, dfill[(row - 1) * dof
//							+ col2 - 1] = 1;
//				}
//			}
//		}
//
//		/* He[1] +  V[1]  ->  He[1]-V[1] */
//		reactants[0] = 1, reactants[1] = N + 1, reactants[2] = 3 * N + 1;
//		for (j = 0; j < 3; j++) {
//			row = reactants[j], col1 = reactants[0], col2 = reactants[1];
//			dfill[(row - 1) * dof + col1 - 1] = 1, dfill[(row - 1) * dof + col2
//					- 1] = 1;
//		}
//
//		/*  He[He]-V[V] + He[he] -> He[He+he]-V[V]  */
//		for (He = 1; He < N; He++) {
//			for (V = 1; V < N + 1; V++) {
//				for (he = 1; he < N - He + 1; he++) {
//					reactants[0] = 3 * N + (He - 1) * N + V, reactants[1] = he, reactants[2] =
//							3 * N + (He + he - 1) * N + V;
//					for (j = 0; j < 3; j++) {
//						row = reactants[j], col1 = reactants[0], col2 =
//								reactants[1];
//						dfill[(row - 1) * dof + col1 - 1] = 1, dfill[(row - 1)
//								* dof + col2 - 1] = 1;
//					}
//				}
//			}
//		}
//		/*  He[He]-V[V] + V[v] -> He[He][V+v] */
//		for (He = 1; He < N + 1; He++) {
//			for (V = 1; V < N; V++) {
//				for (v = 1; v < N - V + 1; v++) {
//					reactants[0] = 3 * N + (He - 1) * N + V, reactants[1] = N
//							+ v, reactants[2] = 3 * N + (He - 1) * N + V + v;
//					for (j = 0; j < 3; j++) {
//						row = reactants[j], col1 = reactants[0], col2 =
//								reactants[1];
//						dfill[(row - 1) * dof + col1 - 1] = 1, dfill[(row - 1)
//								* dof + col2 - 1] = 1;
//					}
//				}
//			}
//		}
//
//		/*  He[He]-V[V]  + He[he]-V[v] -> He[He+he][V+v]  */
//		/*  Currently the reaction rates for this are zero */
//		for (He = 1; He < N; He++) {
//			for (V = 1; V < N; V++) {
//				for (he = 1; he < N - He + 1; he++) {
//					for (v = 1; v < N - V + 1; v++) {
//						reactants[0] = 3 * N + (He - 1) * N + V, reactants[1] =
//								3 * N + (he - 1) * N + V, reactants[2] = 3 * N
//								+ (He + he - 1) * N + V + v;
//						for (j = 0; j < 3; j++) {
//							row = reactants[j], col1 = reactants[0], col2 =
//									reactants[1];
//							dfill[(row - 1) * dof + col1 - 1] = 1, dfill[(row
//									- 1) * dof + col2 - 1] = 1;
//						}
//					}
//				}
//			}
//		}
//		/*  V[V] + I[I]  ->   V[V-I] if V > I else I[I-V] */
//		/*  What should the correct reaction rate should be? */
//		for (V = 1; V < N + 1; V++) {
//			for (I = 1; I < V; I++) {
//				reactants[0] = N + V, reactants[1] = 2 * N + I, reactants[2] = N
//						+ V - I;
//				for (j = 0; j < 3; j++) {
//					row = reactants[j], col1 = reactants[0], col2 = reactants[1];
//					dfill[(row - 1) * dof + col1 - 1] = 1, dfill[(row - 1) * dof
//							+ col2 - 1] = 1;
//				}
//			}
//			for (I = V + 1; I < N + 1; I++) {
//				reactants[0] = N + V, reactants[1] = 2 * N + I, reactants[2] = 2
//						* N + I - V;
//				for (j = 0; j < 3; j++) {
//					row = reactants[j], col1 = reactants[0], col2 = reactants[1];
//					dfill[(row - 1) * dof + col1 - 1] = 1, dfill[(row - 1) * dof
//							+ col2 - 1] = 1;
//				}
//			}
//		}
//	}
//	/* -------------------------------------------------------------------------
//	 ---- Compute dissociation terms that removes an item from a cluster
//	 I assume dissociation means losing only a single item from a cluster
//	 I cannot tell from the notes if clusters can break up into any sub-size.
//	 */
//	if (!ctx->nodissociations) {
//		/*   He[He] ->  He[He-1] + He[1] */
//		for (He = 2; He < N + 1; He++) {
//			reactants[0] = He, reactants[1] = He - 1, reactants[2] = 1;
//
//			for (j = 0; j < 3; j++) {
//				row = reactants[j], col1 = reactants[0];
//				dfill[(row - 1) * dof + col1 - 1] = 1;
//			}
//		}
//		/*   V[V] ->  V[V-1] + V[1] */
//		for (V = 2; V < N + 1; V++) {
//			reactants[0] = N + V, reactants[1] = N + V - 1, reactants[2] = N + 1;
//
//			for (j = 0; j < 3; j++) {
//				row = reactants[j], col1 = reactants[0];
//				dfill[(row - 1) * dof + col1 - 1] = 1;
//			}
//		}
//
//		/*   I[I] ->  I[I-1] + I[1] */
//		for (I = 2; I < N + 1; I++) {
//			reactants[0] = 2 * N + I, reactants[1] = 2 * N + I - 1, reactants[2] =
//					2 * N + 1;
//
//			for (j = 0; j < 3; j++) {
//				row = reactants[j], col1 = reactants[0];
//				dfill[(row - 1) * dof + col1 - 1] = 1;
//			}
//		}
//
//		/* He[1]-V[1]  ->  He[1] + V[1] */
//		reactants[0] = 3 * N + 1, reactants[1] = 1, reactants[2] = N + 1;
//
//		for (j = 0; j < 3; j++) {
//			row = reactants[j], col1 = reactants[0];
//			dfill[(row - 1) * dof + col1 - 1] = 1;
//		}
//
//		/*   He[He]-V[1] ->  He[He] + V[1]  */
//		for (He = 2; He < N + 1; He++) {
//			reactants[0] = 3 * N + (He - 1) * N + 1, reactants[1] = He, reactants[2] =
//					N + 1;
//
//			for (j = 0; j < 3; j++) {
//				row = reactants[j], col1 = reactants[0];
//				dfill[(row - 1) * dof + col1 - 1] = 1;
//			}
//		}
//
//		/*   He[1]-V[V] ->  He[1] + V[V]  */
//		for (V = 2; V < N + 1; V++) {
//			reactants[0] = 3 * N + V, reactants[1] = 1, reactants[2] = N + V;
//
//			for (j = 0; j < 3; j++) {
//				row = reactants[j], col1 = reactants[0];
//				dfill[(row - 1) * dof + col1 - 1] = 1;
//			}
//		}
//
//		/*   He[He]-V[V] ->  He[He-1]-V[V] + He[1]  */
//		for (He = 2; He < N + 1; He++) {
//			for (V = 2; V < N + 1; V++) {
//				reactants[0] = 3 * N + (He - 1) * N + V, reactants[1] = 3 * N
//						+ (He - 2) * N + V, reactants[2] = 1;
//
//				for (j = 0; j < 3; j++) {
//					row = reactants[j], col1 = reactants[0];
//					dfill[(row - 1) * dof + col1 - 1] = 1;
//				}
//			}
//		}
//
//		/*   He[He]-V[V] ->  He[He]-V[V-1] + V[1]  */
//		for (He = 2; He < N + 1; He++) {
//			for (V = 2; V < N + 1; V++) {
//				reactants[0] = 3 * N + (He - 1) * N + V, reactants[1] = 3 * N
//						+ (He - 1) * N + V - 1, reactants[2] = N + 1;
//
//				for (j = 0; j < 3; j++) {
//					row = reactants[j], col1 = reactants[0];
//					dfill[(row - 1) * dof + col1 - 1] = 1;
//				}
//			}
//		}
//
//		/*   He[He]-V[V] ->  He[He]-V[V+1] + I[1]  */
//		/* Again, what is the reasonable dissociation rate? */
//		for (He = 1; He < N + 1; He++) {
//			for (V = 1; V < N; V++) {
//				reactants[0] = 3 * N + (He - 1) * N + V, reactants[1] = 3 * N
//						+ (He - 1) * N + V + 1, reactants[2] = 2 * N + 1;
//
//				for (j = 0; j < 3; j++) {
//					row = reactants[j], col1 = reactants[0];
//					dfill[(row - 1) * dof + col1 - 1] = 1;
//				}
//			}
//		}
//	}

	return 0;
}

//! The Constructor
PetscSolver::PetscSolver() {
	numCLIArgs = 0;
	CLIArgs = NULL;
}

//! The Destructor
PetscSolver::~PetscSolver() {
}

/**
 * This operation transfers the input arguments passed to the program on
 * startup to the solver. These options are static options specified at
 * the start of the program whereas the options passed to setOptions() may
 * change.
 * @param argc The number of command line arguments
 * @param argv The array of command line arguments
 */
void PetscSolver::setCommandLineOptions(int argc, char **argv) {

	numCLIArgs = argc;
	CLIArgs = argv;
}

/**
 * This operation sets the PSIClusterNetworkLoader that should be used by
 * the ISolver to load the ReactionNetwork.
 * @param networkLoader The PSIClusterNetworkLoader that will load the
 * network.
 */
void PetscSolver::setNetworkLoader(
		std::shared_ptr<PSIClusterNetworkLoader> networkLoader) {

// Store the loader and load the network
	this->networkLoader = networkLoader;
	network = networkLoader->load();

// Debug
	std::cout << "PETScSolver Message: " << "Loaded network of size "
			<< network->size() << "." << std::endl;

	return;
}

/**
 * This operation sets the run-time options of the solver. The map is a set
 * of key-value std::string pairs that are interpreted by the solver. These
 * options may change during execution, but it is up to Solvers to monitor
 * the map for changes and they may do so at their discretion.
 * @param options The set of options as key-value pairs with option names
 * for keys and associated values mapped to those keys. A relevant example
 * is "startTime" and "0.01" where both are of type std::string.
 */
void PetscSolver::setOptions(std::map<std::string, std::string> options) {
}

/**
 * This operation sets up the mesh that will be used by the solver and
 * initializes the data on that mesh. This operation will throw an exception
 * of type std::string if the mesh can not be setup.
 */
void PetscSolver::setupMesh() {
}

#undef __FUNCT__
#define __FUNCT__ "initialize"

/**
 * This operation performs all necessary initialization for the solver
 * possibly including but not limited to setting up MPI and loading initial
 * conditions. If the solver can not be initialized, this operation will
 * throw an exception of type std::string.
 */
void PetscSolver::initialize() {

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Initialize program
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	PetscInitialize(&numCLIArgs, &CLIArgs, (char*) 0, help);

	return;
}

/**
 * This operation directs the Solver to perform the solve. If the solve
 * fails, it will throw an exception of type std::string.
 */
void PetscSolver::solve() {

	// Check the network before getting busy.
	if (!network) {
		throw std::string("PetscSolver Exception: Network not set!");
	}

	// The degrees of freedom should be equal to the number of reactants.
	dof = network->size();

	PetscFunctionBeginUser;
	ctx.noreactions = PETSC_FALSE;
	ctx.nodissociations = PETSC_FALSE;

	ierr = PetscOptionsHasName(NULL, "-noreactions", &ctx.noreactions);
	checkPetscError(ierr);
	ierr = PetscOptionsHasName(NULL, "-nodissociations", &ctx.nodissociations);
	checkPetscError(ierr);

	ctx.HeDiffusion[1] = 1000 * 2.95e-4; /* From Tibo's notes times 1,000 */
	ctx.HeDiffusion[2] = 1000 * 3.24e-4;
	ctx.HeDiffusion[3] = 1000 * 2.26e-4;
	ctx.HeDiffusion[4] = 1000 * 1.68e-4;
	ctx.HeDiffusion[5] = 1000 * 5.20e-5;
	ctx.VDiffusion[1] = 1000 * 2.71e-3;
	ctx.IDiffusion[1] = 1000 * 2.13e-4;
	ctx.forcingScale = 100.; /* made up numbers */
	ctx.reactionScale = .001;
	ctx.dissociationScale = .0001;
	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create distributed array (DMDA) to manage parallel grid and vectors
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMDACreate1d(PETSC_COMM_WORLD, DMDA_BOUNDARY_MIRROR, -8, dof, 1,
			NULL, &da);
	checkPetscError(ierr);

	/* The only spatial coupling in the Jacobian (diffusion) is for the first 5 He, the first V, and the first I.
	 The ofill (thought of as a dof by dof 2d (row-oriented) array represents the nonzero coupling between degrees
	 of freedom at one point with degrees of freedom on the adjacent point to the left or right. A 1 at i,j in the
	 ofill array indicates that the degree of freedom i at a point is coupled to degree of freedom j at the
	 adjacent point. In this case ofill has only a few diagonal entries since the only spatial coupling is regular diffusion. */
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &ofill);
	checkPetscError(ierr);
	ierr = PetscMalloc(dof * dof * sizeof(PetscInt), &dfill);
	checkPetscError(ierr);
	ierr = PetscMemzero(ofill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr);
	ierr = PetscMemzero(dfill, dof * dof * sizeof(PetscInt));
	checkPetscError(ierr);

	for (He = 0; He < PetscMin(N, 5); He++)
		ofill[He * dof + He] = 1;
	ofill[N * dof + N] = ofill[2 * N * dof + 2 * N] = 1;

	ierr = DMDASetBlockFills(da, NULL, ofill);
	checkPetscError(ierr);
	ierr = PetscFree(ofill);
	checkPetscError(ierr);
	ierr = getDiagonalFill(dfill, dof * dof);
	checkPetscError(ierr);
	ierr = PetscFree(dfill);
	checkPetscError(ierr);

	/*  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Extract global vector from DMDA to hold solution
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = DMCreateGlobalVector(da, &C);
	checkPetscError(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Create timestepping solver context
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = TSCreate(PETSC_COMM_WORLD, &ts);
	checkPetscError(ierr);
	ierr = TSSetType(ts, TSARKIMEX);
	checkPetscError(ierr);
	ierr = TSSetDM(ts, da);
	checkPetscError(ierr);
	ierr = TSSetProblemType(ts, TS_NONLINEAR);
	checkPetscError(ierr);
	ierr = TSSetIFunction(ts, NULL, IFunction, &ctx);
	checkPetscError(ierr);
	ierr = TSSetSolution(ts, C);
	checkPetscError(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set solver options
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = TSSetInitialTimeStep(ts, 0.0, .001);
	checkPetscError(ierr);
	ierr = TSSetDuration(ts, 100, 50.0);
	checkPetscError(ierr);
	ierr = TSSetFromOptions(ts);
	checkPetscError(ierr);
	ierr = setupPetscMonitor(ts);
	checkPetscError(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Set initial conditions
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = setupInitialConditions(da, C);
	checkPetscError(ierr);

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Solve the ODE system
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	if (ts != NULL && C != NULL) {
		ierr = TSSolve(ts, C);
		checkPetscError(ierr);
	} else {
		throw std::string(
				"PetscSolver Exception: Unable to solve! Data not configured properly.");
	}

}

/**
 * This operation performs all necessary finalization for the solver
 * including but not limited to cleaning up memory, finalizing MPI and
 * printing diagnostic information. If the solver can not be finalized,
 * this operation will throw an exception of type std::string.
 */
void PetscSolver::finalize() {

	/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	 Free work space.
	 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
	ierr = VecDestroy(&C);
	checkPetscError(ierr);
	ierr = TSDestroy(&ts);
	checkPetscError(ierr);
	ierr = DMDestroy(&da);
	checkPetscError(ierr);
	ierr = PetscFinalize();
	if (petscReturn() != 0) {
		throw std::string("PetscSolver Exception: Unable to finalize solve!");
	}

}

} /* end namespace xolotlSolver */
