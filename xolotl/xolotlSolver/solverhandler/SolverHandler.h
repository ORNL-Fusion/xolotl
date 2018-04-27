#ifndef SOLVERHANDLER_H
#define SOLVERHANDLER_H

// Includes
#include "ISolverHandler.h"
#include <HDF5Utils.h>
#include "RandomNumberGenerator.h"

namespace xolotlSolver {

/**
 * This class and its subclasses realize the ISolverHandler interface to solve the
 * advection-diffusion-reaction problem with currently supported solvers.
 */
class SolverHandler: public ISolverHandler {
protected:

	//! The name of the network file
	std::string networkName;

	//! The original network created from the network loader.
	xolotlCore::IReactionNetwork& network;

	//! Vector storing the grid in the x direction
	std::vector<double> grid;

	//! The number of grid points in the depth direction.
	int nX;

	//! The number of grid points in the Y direction.
	int nY;

	//! The number of grid points in the Z direction.
	int nZ;

	//! The grid step size in the depth direction.
	double hX;

	//! The grid step size in the y direction.
	double hY;

	//! The grid step size in the z direction.
	double hZ;

	//! The number of grid points by which the boundary condition should be shifted at the left side.
	int leftOffset;

	//! The number of grid points by which the boundary condition should be shifted at the right side.
	int rightOffset;

	//! The initial vacancy concentration.
	double initialVConc;

	//! The original flux handler created.
	xolotlCore::IFluxHandler *fluxHandler;

	//! The original temperature handler created.
	xolotlCore::ITemperatureHandler *temperatureHandler;

	//! The original diffusion handler created.
	xolotlCore::IDiffusionHandler *diffusionHandler;

	//! The vector of advection handlers.
	std::vector<xolotlCore::IAdvectionHandler *> advectionHandlers;

	//! The original modified trap-mutation handler created.
	xolotlCore::ITrapMutationHandler *mutationHandler;

	//! The number of dimensions for the problem.
	int dimension;

	//! The portion of void at the beginning of the problem.
	double portion;

	//! If the user wants to use a regular grid.
	bool useRegularGrid;

	//! If the user wants to move the surface.
	bool movingSurface;

	//! If the user wants to burst bubbles.
	bool bubbleBursting;

	//! The sputtering yield for the problem.
	double sputteringYield;

	//! The depth parameter for the bubble bursting.
	double tauBursting;

    //! The value to use to seed the random number generator.
    unsigned int rngSeed;

    //! The random number generator to use.
    std::unique_ptr<RandomNumberGenerator<int, unsigned int>> rng;

	//! Method generating the grid in the x direction
	void generateGrid(int nx, double hx, int surfacePos) {
		// Clear the grid
		grid.clear();

		// Check if the user wants a regular grid
		if (useRegularGrid) {
			// The grid will me made of nx + 1 points separated by hx nm
			for (int l = 0; l <= nx + 1; l++) {
				grid.push_back((double) l * hx);
			}
		}
		// If it is not regular do a fine mesh close to the surface and
		// increase the step size when away from the surface
		else {
			// Initialize the value of the previous point
			double previousPoint = 0.0;

			// Loop on all the grid points
			for (int l = 0; l <= nx + 1; l++) {
				// Add the previous point
				grid.push_back(previousPoint);
				// 0.1nm step near the surface (x < 2.5nm)
				if (l < surfacePos + 26) {
					previousPoint += 0.1;
				}
				// Then 0.25nm (2.5nm < x < 5.0nm)
				else if (l < surfacePos + 36) {
					previousPoint += 0.25;
				}
				// Then 0.5nm (5.0nm < x < 7.5nm)
				else if (l < surfacePos + 41) {
					previousPoint += 0.5;
				}
				// Then 1.0nm step size (7.5nm < x < 50.5)
				else if (l < surfacePos + 84) {
					previousPoint += 1.0;
				}
				// Then 2.0nm step size (50.5nm < x < 100.5)
				else if (l < surfacePos + 109) {
					previousPoint += 2.0;
				}
				// Then 5.0nm step size (100.5nm < x < 150.5)
				else if (l < surfacePos + 119) {
					previousPoint += 5.0;
				}
				// Then 10.0nm step size (150.5nm < x < 300.5)
				else if (l < surfacePos + 134) {
					previousPoint += 10.0;
				}
				// Then 20.0nm step size (300.5nm < x < 500.5)
				else if (l < surfacePos + 144) {
					previousPoint += 20.0;
				}
				// Then 50.0nm step size (500.5nm < x < 1000.5)
				else if (l < surfacePos + 154) {
					previousPoint += 50.0;
				}
				// Then 100.0nm step size (1000.5nm < x < 5000.5)
				else if (l < surfacePos + 194) {
					previousPoint += 100.0;
				}
				// Then 200.0nm step size (5000.5nm < x < 10000.5)
				else if (l < surfacePos + 219) {
					previousPoint += 200.0;
				}
				// Then 500.0nm step size (10000.5nm < x < 20000.5)
				else if (l < surfacePos + 239) {
					previousPoint += 500.0;
				}
				// Then 1.0um step size (20000.5nm < x)
				else {
					previousPoint += 1000.0;
				}
			}
		}

		return;
	}

	/**
	 * Constructor.
	 *
	 * @param _network The reaction network to use.
	 */
	SolverHandler(xolotlCore::IReactionNetwork& _network) :
			network(_network), networkName(""), nX(0), nY(0), nZ(0), hX(0.0), hY(
					0.0), hZ(0.0), leftOffset(0), rightOffset(0), initialVConc(
					0.0), dimension(-1), portion(0.0), useRegularGrid(true), movingSurface(
					false), bubbleBursting(false), sputteringYield(0.0), fluxHandler(
					nullptr), temperatureHandler(nullptr), diffusionHandler(
					nullptr), mutationHandler(nullptr), tauBursting(10.0) {
	}

public:

	//! The Constructor
	SolverHandler() = delete;

	~SolverHandler() {
	}

	/**
	 * Initialize all the physics handlers that are needed to solve the ADR equations.
	 * \see ISolverHandler.h
	 */
	void initializeHandlers(
			std::shared_ptr<xolotlFactory::IMaterialFactory> material,
			std::shared_ptr<xolotlCore::ITemperatureHandler> tempHandler,
			const xolotlCore::Options &options) override {

        // Determine who I am.
        int myProcId = -1;
        MPI_Comm_rank(MPI_COMM_WORLD, &myProcId);

        // Initialize our random number generator.
        bool useRNGSeedFromOptions = false;
        bool printRNGSeed = false;
        std::tie(useRNGSeedFromOptions, rngSeed) = options.getRNGSeed();
        if(not useRNGSeedFromOptions) {
            // User didn't give a seed value to use, so
            // use something based on current time and our proc id
            // so that it is different from run to run, and should
            // be different across all processes within a given run.
            rngSeed = time(NULL) + myProcId;
        }
        if(options.printRNGSeed()) {
            std::cout << "Proc " << myProcId
                << " using RNG seed value " << rngSeed
                << std::endl;
        }
        rng = std::unique_ptr<RandomNumberGenerator<int, unsigned int>>(new RandomNumberGenerator<int, unsigned int>(rngSeed));

		// Set the network loader
		networkName = options.getNetworkFilename();

		// Set the grid options
		if (options.useHDF5()) {
			// Get starting conditions from HDF5 file
			int nx = 0, ny = 0, nz = 0;
			double hx = 0.0, hy = 0.0, hz = 0.0;
			xolotlCore::HDF5Utils::readHeader(networkName, nx, hx, ny, hy, nz,
					hz);
			nX = nx, nY = ny, nZ = nz;
			hX = hx, hY = hy, hZ = hz;
		} else {
			nX = options.getNX(), nY = options.getNY(), nZ = options.getNZ();
			hX = options.getXStepSize(), hY = options.getYStepSize(), hZ =
					options.getZStepSize();
		}

		// Set the flux handler
		fluxHandler =
				(xolotlCore::IFluxHandler *) material->getFluxHandler().get();

		// Set the temperature handler
		temperatureHandler =
				(xolotlCore::ITemperatureHandler *) tempHandler.get();

		// Set the diffusion handler
		diffusionHandler =
				(xolotlCore::IDiffusionHandler *) material->getDiffusionHandler().get();

		// Set the advection handlers
		auto handlers = material->getAdvectionHandler();
		for (auto handler : handlers) {
			advectionHandlers.push_back(handler.get());
		}

		// Set the modified trap-mutation handler
		mutationHandler =
				(xolotlCore::ITrapMutationHandler *) material->getTrapMutationHandler().get();

		// Set the initial vacancy concentration
		initialVConc = options.getInitialVConcentration();

		// Set the number of dimension
		dimension = options.getDimensionNumber();

		// Set the void portion
		portion = options.getVoidPortion();

		// Set the sputtering yield
		sputteringYield = options.getSputteringYield();

		// Set the sputtering yield
		tauBursting = options.getBurstingDepth();

		// Look at if the user wants to use a regular grid in the x direction
		useRegularGrid = options.useRegularXGrid();

		// Set the boundary conditions (= 1: free surface; = 0: mirror)
		leftOffset = options.getLeftBoundary();
		rightOffset = options.getRightBoundary();

		// Should we be able to move the surface?
		auto map = options.getProcesses();
		movingSurface = map["movingSurface"];
		// Should we be able to burst bubble?
		bubbleBursting = map["bursting"];

		// Some safeguards about what to use with what
		if (leftOffset == 0
				&& (map["advec"] || map["modifiedTM"] || map["movingSurface"]
						|| map["bursting"])) {
			throw std::string(
					"\nThe left side of the grid is set to use a reflective boundary condition "
							"but you want to use processes that are intrinsically related to "
							"a free surface (advection, modified trap mutation, moving surface, bubble bursting).");
		}

		// Complains if processes that should not be used together are used
		if (map["attenuation"] && !map["modifiedTM"]) {
			throw std::string(
					"\nYou want to use the attenuation on the modified trap mutation "
							"but you are not using the modifiedTM process, it doesn't make any sense.");
		}
		if (map["modifiedTM"] && !map["reaction"]) {
			throw std::string(
					"\nYou want to use the modified trap mutation but the reaction process is not set,"
							" it doesn't make any sense.");
		}

		return;
	}

	/**
	 * Get the grid in the x direction.
	 * \see ISolverHandler.h
	 */
	std::vector<double> getXGrid() const override {
		return grid;
	}

	/**
	 * Get the step size in the y direction.
	 * \see ISolverHandler.h
	 */
	double getStepSizeY() const override {
		return hY;
	}

	/**
	 * Get the step size in the z direction.
	 * \see ISolverHandler.h
	 */
	double getStepSizeZ() const override {
		return hZ;
	}

	/**
	 * Get the number of dimensions of the problem.
	 * \see ISolverHandler.h
	 */
	int getDimension() const override {
		return dimension;
	}

	/**
	 * Get the initial vacancy concentration.
	 * \see ISolverHandler.h
	 */
	double getInitialVConc() const override {
		return initialVConc;
	}

	/**
	 * Get the sputtering yield.
	 * \see ISolverHandler.h
	 */
	double getSputteringYield() const override {
		return sputteringYield;
	}

	/**
	 * Get the depth parameter for bursting.
	 * \see ISolverHandler.h
	 */
	double getTauBursting() const override {
		return tauBursting;
	}

	/**
	 * Get the grid right offset.
	 * \see ISolverHandler.h
	 */
	int getRightOffset() const {
		return rightOffset;
	}

	/**
	 * To know if the surface should be able to move.
	 * \see ISolverHandler.h
	 */
	bool moveSurface() const override {
		return movingSurface;
	}

	/**
	 * To know if the bubble bursting should be used.
	 * \see ISolverHandler.h
	 */
	bool burstBubbles() const override {
		return bubbleBursting;
	}

	/**
	 * Get the flux handler.
	 * \see ISolverHandler.h
	 */
	xolotlCore::IFluxHandler *getFluxHandler() const override {
		return fluxHandler;
	}

	/**
	 * Get the temperature handler.
	 * \see ISolverHandler.h
	 */
	xolotlCore::ITemperatureHandler *getTemperatureHandler() const override {
		return temperatureHandler;
	}

	/**
	 * Get the advection handler.
	 * \see ISolverHandler.h
	 */
	xolotlCore::IAdvectionHandler *getAdvectionHandler() const override {
		return advectionHandlers[0];
	}

	/**
	 * Get the advection handlers.
	 * \see ISolverHandler.h
	 */
	std::vector<xolotlCore::IAdvectionHandler *> getAdvectionHandlers() const override {
		return advectionHandlers;
	}

	/**
	 * Get the modified trap-mutation handler.
	 * \see ISolverHandler.h
	 */
	xolotlCore::ITrapMutationHandler *getMutationHandler() const override {
		return mutationHandler;
	}

	/**
	 * Get the network.
	 * \see ISolverHandler.h
	 */
	xolotlCore::IReactionNetwork& getNetwork() const override {
		return network;
	}

	/**
	 * Get the network name.
	 * \see ISolverHandler.h
	 */
	std::string getNetworkName() const override {
		return networkName;
	}

    /**
     * Access the random number generator
     * The generator will have already been seeded.
     *
     * @return The RandomNumberGenerator object to use.
     */
    RandomNumberGenerator<int, unsigned int>& getRNG(void) const override {
        return *rng;
    }
}
;
//end class SolverHandler

} /* end namespace xolotlSolver */
#endif
