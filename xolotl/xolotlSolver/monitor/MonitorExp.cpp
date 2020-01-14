// Includes
#include "PetscSolver.h"
#include <xolotlPerf.h>
#include <Constants.h>
#include <petscts.h>
#include <petscsys.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <memory>
#include "xolotlSolver/monitor/Monitor.h"

namespace xolotlSolver {

// Declaration of the functions defined in Monitor.cpp
extern PetscErrorCode monitorTime(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);
extern PetscErrorCode computeFluence(TS ts, PetscInt timestep, PetscReal time,
		Vec solution, void *ictx);

// Declaration of the variables defined in Monitor.cpp
extern double previousTime;

// Declare the vector that will store the Id of the helium clusters
std::vector<int> indicesExp, heIndicesExp, dIndicesExp, tIndicesExp,
		vIndicesExp, iIndicesExp;
// Declare the vector that will store the weight of the helium clusters
// (their He composition)
std::vector<int> weightsExp, heWeightsExp, dWeightsExp, tWeightsExp,
		vWeightsExp, iWeightsExp;
// Declare the vector that will store the radii of the helium clusters
std::vector<double> radiiExp;

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeXenonRetentionExp")
/**
 * This is a monitoring method that will compute the xenon retention
 */
PetscErrorCode computeXenonRetentionExp(TS ts, PetscInt, PetscReal time,
		Vec solution, void *) {
	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration and other values over the grid
	double xeConcentration = 0.0, bubbleConcentration = 0.0, radii = 0.0,
			partialBubbleConcentration = 0.0, partialRadii = 0.0, partialSize =
					0.0;

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Get the minimum size for the radius
	auto minSizes = solverHandler.getMinSizes();

	// Loop on all the indices
	for (unsigned int i = 0; i < indicesExp.size(); i++) {
		// Add the current concentration times the number of xenon in the cluster
		// (from the weight vector)
		double conc = gridPointSolution[indicesExp[i]];
		xeConcentration += conc * weightsExp[i];
		bubbleConcentration += conc;
		radii += conc * radiiExp[i];
		if (weightsExp[i] >= minSizes[0] && conc > 1.0e-16) {
			partialBubbleConcentration += conc;
			partialRadii += conc * radiiExp[i];
			partialSize += conc * weightsExp[i];
		}
	}

	// Print the result
	std::cout << "\nTime: " << time << std::endl;
	std::cout << "Xenon concentration = " << xeConcentration << std::endl
			<< std::endl;

	// Make sure the average partial radius makes sense
	double averagePartialRadius = partialRadii / partialBubbleConcentration;
//	double minRadius = pow(
//			(3.0 * (double) minSizes[0])
//					/ (4.0 * xolotlCore::pi * network.getDensity()),
//			(1.0 / 3.0));
//	if (partialBubbleConcentration < 1.e-16 || averagePartialRadius < minRadius)
//		averagePartialRadius = minRadius;

	// Uncomment to write the retention and the fluence in a file
	std::ofstream outputFile;
	outputFile.open("retentionOut.txt", ios::app);
	outputFile << time << " " << xeConcentration << " "
			<< radii / bubbleConcentration << " " << averagePartialRadius << " "
			<< partialBubbleConcentration << " "
			<< partialSize / partialBubbleConcentration << std::endl;
	outputFile.close();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ Actual__FUNCT__("xolotlSolver", "computeHeliumRetentionExp")
/**
 * This is a monitoring method that will compute the helium retention
 */
PetscErrorCode computeHeliumRetentionExp(TS ts, PetscInt, PetscReal time,
		Vec solution, void *) {
	// Initial declarations
	PetscErrorCode ierr;

	PetscFunctionBeginUser;

	// Get the da from ts
	DM da;
	ierr = TSGetDM(ts, &da);
	CHKERRQ(ierr);

	// Get the array of concentration
	PetscReal **solutionArray;
	ierr = DMDAVecGetArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	// Store the concentration and other values over the grid
	double heConc = 0.0, dConc = 0.0, tConc = 0.0, vConc = 0.0, iConc = 0.0;

	// Declare the pointer for the concentrations at a specific grid point
	PetscReal *gridPointSolution;

	// Get the pointer to the beginning of the solution data for this grid point
	gridPointSolution = solutionArray[0];

	// Loop on all the indices
	for (unsigned int i = 0; i < heIndicesExp.size(); i++) {
		// Add the current concentration times the number of helium in the cluster
		// (from the weight vector)
		heConc += gridPointSolution[heIndicesExp[i]] * heWeightsExp[i];
	}
	// D
	for (unsigned int i = 0; i < dIndicesExp.size(); i++) {
		dConc += gridPointSolution[dIndicesExp[i]] * dWeightsExp[i];
	}
	// T
	for (unsigned int i = 0; i < tIndicesExp.size(); i++) {
		tConc += gridPointSolution[tIndicesExp[i]] * tWeightsExp[i];
	}
	// V
	for (unsigned int i = 0; i < vIndicesExp.size(); i++) {
		vConc += gridPointSolution[vIndicesExp[i]] * vWeightsExp[i];
	}
	// I
	for (unsigned int i = 0; i < iIndicesExp.size(); i++) {
		iConc += gridPointSolution[iIndicesExp[i]] * iWeightsExp[i];
	}

	// Print the result
	std::cout << "\nTime: " << time << std::endl;
	std::cout << "Helium concentration = " << heConc << std::endl;
	std::cout << "Deuterium concentration = " << dConc << std::endl;
	std::cout << "Tritium concentration = " << tConc << std::endl;
	std::cout << "Vacancy concentration = " << vConc << std::endl;
	std::cout << "Interstital concentration = " << iConc << std::endl
			<< std::endl;

	// Uncomment to write the retention and the fluence in a file
	std::ofstream outputFile;
	outputFile.open("retentionOut.txt", ios::app);
	outputFile << time << " " << heConc << " " << dConc << " " << tConc << " "
			<< vConc << " " << iConc << std::endl;
	outputFile.close();

	// Restore the solutionArray
	ierr = DMDAVecRestoreArrayDOFRead(da, solution, &solutionArray);
	CHKERRQ(ierr);

	PetscFunctionReturn(0);
}

/**
 * This operation sets up different monitors
 *  depending on the options.
 * @param ts The time stepper
 * @return A standard PETSc error code
 */
PetscErrorCode setupPetscExpMonitor(TS ts) {
	PetscErrorCode ierr;

	// Flags to launch the monitors or not
	PetscBool flagXeRetention, flagHeRetention;

	// Check the option -xenon_retention
	ierr = PetscOptionsHasName(NULL, NULL, "-xenon_retention",
			&flagXeRetention);
	checkPetscError(ierr,
			"setupPetscExpMonitor: PetscOptionsHasName (-xenon_retention) failed.");

	// Check the option -helium_retention
	ierr = PetscOptionsHasName(NULL, NULL, "-helium_retention",
			&flagHeRetention);
	checkPetscError(ierr,
			"setupPetscExpMonitor: PetscOptionsHasName (-helium_retention) failed.");

	// Get the solver handler
	auto& solverHandler = PetscSolver::getSolverHandler();

	// Get the network and its size
	auto& network = solverHandler.getExpNetwork();
	const int networkSize = network.getDOF();

	// Set the monitor to compute the xenon fluence and the retention
	// for the retention calculation
	if (flagXeRetention) {

		using NetworkType =
		experimental::NEReactionNetwork;
		using Spec = typename NetworkType::Species;
		using Composition = typename NetworkType::Composition;

		for (std::size_t i = 0; i < networkSize; ++i) {
			const auto& cl1Reg = network.getCluster(i).getRegion();
			Composition lo1 = cl1Reg.getOrigin();
			// Add the Id to the vector
			indicesExp.push_back(i);
			// Add the number of xenon of this cluster to the weight
			weightsExp.push_back(lo1[Spec::Xe]);
			radiiExp.push_back(network.getCluster(i).getReactionRadius());
		}

		// computeXenonRetentionExp will be called at each timestep
		ierr = TSMonitorSet(ts, computeXenonRetentionExp, NULL, NULL);
		checkPetscError(ierr,
				"setupPetscExpMonitor: TSMonitorSet (computeXenonRetentionExp) failed.");

		// Uncomment to clear the file where the retention will be written
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt");
		outputFile.close();
	}

	// Set the monitor to compute the helium fluence and the retention
	// for the retention calculation
	if (flagHeRetention) {
//		using NetworkType =
//		experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
//		using Spec = typename NetworkType::Species;
//		using Composition = typename NetworkType::Composition;
//
//		// Loop on the clusters
//		for (std::size_t i = 0; i < networkSize; ++i) {
//			const auto& cl1Reg = network.getCluster(i).getRegion();
//			Composition lo1 = cl1Reg.getOrigin();
//			if (lo1[Spec::He] > 0) {
//				heIndicesExp.push_back(i);
//				heWeightsExp.push_back(lo1[Spec::He]);
//			}
//			if (lo1[Spec::D] > 0) {
//				dIndicesExp.push_back(i);
//				dWeightsExp.push_back(lo1[Spec::D]);
//			}
//			if (lo1[Spec::T] > 0) {
//				tIndicesExp.push_back(i);
//				tWeightsExp.push_back(lo1[Spec::T]);
//			}
//			if (lo1[Spec::V] > 0) {
//				vIndicesExp.push_back(i);
//				vWeightsExp.push_back(lo1[Spec::V]);
//			}
//			if (lo1[Spec::I] > 0) {
//				iIndicesExp.push_back(i);
//				iWeightsExp.push_back(lo1[Spec::I]);
//			}
//		}

		// computeXenonRetentionExp will be called at each timestep
		ierr = TSMonitorSet(ts, computeHeliumRetentionExp, NULL, NULL);
		checkPetscError(ierr,
				"setupPetscExpMonitor: TSMonitorSet (computeHeliumRetentionExp) failed.");

		// Uncomment to clear the file where the retention will be written
		std::ofstream outputFile;
		outputFile.open("retentionOut.txt");
		outputFile.close();
	}

	// Set the monitor to simply change the previous time to the new time
	// monitorTime will be called at each timestep
	ierr = TSMonitorSet(ts, monitorTime, NULL, NULL);
	checkPetscError(ierr,
			"setupPetscExpMonitor: TSMonitorSet (monitorTime) failed.");

	PetscFunctionReturn(0);
}

}

/* end namespace xolotlSolver */
