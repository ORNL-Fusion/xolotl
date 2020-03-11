// Includes
#include "XolotlDLinterface.h"
#include <iostream>
#include <Options.h>
#include <IMaterialFactory.h>
#include <TemperatureHandlerFactory.h>
#include <VizHandlerRegistryFactory.h>
#include <xolotlPerf.h>
#include <SolverHandlerFactory.h>
#include <IReactionHandlerFactory.h>
#include <ctime>
#include <MPIUtils.h>

class XolotlInterface : public XolotlDLinterface {
public:
	void printSomething() {
		std::cout << "I'm in Xolotl !!!" << std::endl;
		return;
	}

	//std::shared_ptr<xolotlSolver::PetscSolver> initializeXolotl(
	//	int argc, char **argv, MPI_Comm comm) {
  std::shared_ptr<xolotlSolver::PetscSolver> initializeXolotl(
		int argc, char **argv, MPI_Comm comm, bool isStandalone) {
		// Local Declarations
		std::shared_ptr<xolotlSolver::PetscSolver> solver;

		// Initialize the MPI communicator to use
		xolotlCore::MPIUtils::initialize(comm);
		auto xolotlComm = xolotlCore::MPIUtils::getMPIComm();

		// Get the MPI rank
		int rank;
		MPI_Comm_rank(xolotlComm, &rank);

		if (rank == 0) {
			// Print the start message
			std::cout << "Starting Xolotl Plasma-Surface Interactions Simulator"
				<< std::endl;
			// TODO! Print copyright message
			// Print date and time
			std::time_t currentTime = std::time(NULL);
			std::cout << std::asctime(std::localtime(&currentTime));
		}

		try {
			// Skip the executable name before parsing
			//argc -= 1; // one for the executable name
			//argv += 1; // one for the executable name
			Options opts;
			opts.readParams(argc,argv);
			if (!opts.shouldRun()) {
				std::cerr << "Unable to read the options.  Aborting" << std::endl;
			}

			// Set up our performance data infrastructure.
			xolotlPerf::initialize(opts.getPerfHandlerType());
			auto handlerRegistry = xolotlPerf::getHandlerRegistry();

			//// Get the MPI rank
			//int rank;
			//MPI_Comm_rank(xolotlComm, &rank);

			// Create the material factory
			auto materialFactory =
				xolotlFactory::IMaterialFactory::createMaterialFactory(
						opts.getMaterial(), opts.getDimensionNumber());
			// Initialize it with the options
			materialFactory->initializeMaterial(opts);

			// Initialize the temperature
			if (!xolotlFactory::initializeTempHandler(opts)) {
				std::cerr << "Unable to initialize requested temperature.  Aborting"
					<< std::endl;
			}
			// Get the temperature handler
			auto tempHandler = xolotlFactory::getTemperatureHandler();

			// Initialize the visualization
			if (!xolotlFactory::initializeVizHandler(
				opts.useVizStandardHandlers())) {
					std::cerr
					<< "Unable to initialize requested visualization infrastructure. "
					<< "Aborting" << std::endl;
			}

			// Create the network handler factory
			auto networkFactory =
				xolotlFactory::IReactionHandlerFactory::createNetworkFactory(
						opts.getMaterial());
			// Build a reaction network
			networkFactory->initializeReactionNetwork(opts, handlerRegistry);
			auto& network = networkFactory->getNetworkHandler();

			// Initialize and get the solver handler
			if (!xolotlFactory::initializeDimension(opts, network)) {
				std::cerr << "Unable to initialize dimension from inputs. "
					<< "Aborting" << std::endl;
			}
			auto& solvHandler = xolotlFactory::getSolverHandler();
			// Initialize the solver handler
			solvHandler.initializeHandlers(materialFactory, tempHandler, opts);

			// Setup the solver
			solver = std::shared_ptr<xolotlSolver::PetscSolver>(
				new xolotlSolver::PetscSolver(solvHandler, handlerRegistry));
                        // Initialize the solver
			//solver->setCommandLineOptions(opts.getPetscArgc(), opts.getPetscArgv());
			solver->setCommandLineOptions(opts.getPetscArg());
			//solver->initialize();
			solver->initialize(isStandalone);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return solver;
	}

	void getXolotlGlobalGridInfo(int *dim, bool *regulargrid, int *nx, int *ny, int *nz, double *hx, double *hy, double *hz,
				double *lx, double *ly, double *lz,
				int argc, char **argv){
		std::shared_ptr<xolotlSolver::PetscSolver> solver;
		try {
			// Skip the executable name before parsing
			//argc -= 1; // one for the executable name
			//argv += 1; // one for the executable name

			// Get Option object
			Options opts;
			opts.readParams(argc, argv);
			if (!opts.shouldRun()) {
				std::cerr << "Unable to read the options.  Aborting" << std::endl;
			}

			// Get Xolotl grid dimension
			*dim = opts.getDimensionNumber();

			// Using regular grid?
			*regulargrid = opts.useRegularXGrid();

			// Get Xolotl the total number of nodes along each axis
			*nx = *dim >= 1 ? opts.getNX() : 1;
			*ny = *dim >= 2 ? opts.getNY() : 1;
			*nz = *dim >= 3 ? opts.getNZ() : 1;

			// Get the grid spacing along each axis
			*hx = *dim >= 1 ? opts.getXStepSize() : 0.0;
			*hy = *dim >= 2 ? opts.getYStepSize() : 0.0;
			*hz = *dim >= 3 ? opts.getZStepSize() : 0.0;

			// Get the total domain length of each axis
			*lx = *dim >= 1 ? (double) (*nx - 1) * *hx : 0.0;
			*ly = *dim >= 2 ? (double) (*ny - 1) * *hy : 0.0;
			*lz = *dim >= 3 ? (double) (*nz - 1) * *hz : 0.0;

			// std::cout<< "getXolotlGlobalGridInfo() called"<<std::endl;
			// std::cout<<"_xolotl_dim = "<< *dim <<std::endl;
		  // std::cout<<"_xolotl_lx = "<< *lx <<std::endl;
		  // std::cout<<"_xolotl_ly = "<< *ly <<std::endl;
		  // std::cout<<"_xolotl_lz = "<< *lz <<std::endl;
		  // std::cout<<"_xolotl_nx, _xolotl_ny, _xolotl_nz = "<<*nx<<","<<*ny<<","<<*nz<<std::endl;
		  // std::cout<<"_xolotl_dx, _xolotl_dy, _xolotl_dz = "<<*hx<<","<<*hy<<","<<*hz<<std::endl;


		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return;
	}

	std::vector<double> getXolotlXgrid(std::shared_ptr<xolotlSolver::PetscSolver> solver){
		std::vector<double> toReturn;
		try {
			// Get the solver handler
			auto& solverHandler = solver->getSolverHandler();
			// Get the non-Regular x-grid
			toReturn = solverHandler.getXGrid();
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return toReturn;
	}

	void setTimes(
		std::shared_ptr<xolotlSolver::PetscSolver> solver, double finalTime,
		double dt) {
		try {
			// Set the time in the solver
			solver->setTimes(finalTime, dt);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}
	}

	double getXolotlTimeInterface(std::shared_ptr<xolotlSolver::PetscSolver> solver) {
		double returnval = 0.0;
		try {
			// get the most recent time that Xolotl solver converged
			returnval = solver->getXolotlTime();
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}
		return returnval;
	}

	void solveXolotl(
		std::shared_ptr<xolotlSolver::PetscSolver> solver) {
		try {
			// Launch the PetscSolver
			solver->solve();
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return;
	}

	// std::vector<std::vector<std::vector<std::tuple<double, double, double, double> > > > * getLocalNE(
	// 		std::shared_ptr<xolotlSolver::PetscSolver> solver) {
	// 	std::vector<std::vector<std::vector<std::tuple<double, double, double, double> > > > * toReturn;
	// 	try {
	// 		// Get the solver handler
	// 		auto& solverHandler = solver->getSolverHandler();
	// 		// Get the rate at this location
	// 		//toReturn = solverHandler.getLocalXeRateCpy();
	// 		toReturn = solverHandler.getLocalNECpy();
	// 	} catch (const std::exception& e) {
	// 		std::cerr << e.what() << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (const std::string& error) {
	// 		std::cerr << error << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (...) {
	// 		std::cerr << "Unrecognized exception seen." << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	}
	//
	// 	return toReturn;
	// }

	// std::vector<std::vector<std::vector<double> > > * getLocalXeRate(
	// 		std::shared_ptr<xolotlSolver::PetscSolver> solver) {
	// 	std::vector<std::vector<std::vector<double> > > * toReturn;
	// 	try {
	// 		// Get the solver handler
	// 		auto& solverHandler = solver->getSolverHandler();
	// 		// Get the rate at this location
	// 		//toReturn = solverHandler.getLocalXeRateCpy();
	// 		toReturn = solverHandler.getLocalXeRate();
	// 	} catch (const std::exception& e) {
	// 		std::cerr << e.what() << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (const std::string& error) {
	// 		std::cerr << error << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (...) {
	// 		std::cerr << "Unrecognized exception seen." << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	}
	//
	// 	return toReturn;
	// }



	double getLocalXeRatePoint(
			std::shared_ptr<xolotlSolver::PetscSolver> solver,
			int i, int j, int k) {
		double toReturn;
		try {
			// Get the solver handler
			auto& solverHandler = solver->getSolverHandler();
			// Get the rate at this location
			toReturn = solverHandler.getXeRatePoint(i, j, k);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return toReturn;
	}

	double getLocalXeFluxPoint(
			std::shared_ptr<xolotlSolver::PetscSolver> solver,
			int i, int j, int k) {
		double toReturn;
		try {
			// Get the solver handler
			auto& solverHandler = solver->getSolverHandler();
			// Get the rate at this location
			toReturn = solverHandler.getXeFluxPoint(i, j, k);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return toReturn;
	}

	double getLocalXeConcPoint(
			std::shared_ptr<xolotlSolver::PetscSolver> solver,
			int i, int j, int k) {
		double toReturn;
		try {
			// Get the solver handler
			auto& solverHandler = solver->getSolverHandler();
			// Get the rate at this location
			toReturn = solverHandler.getXeConcPoint(i, j, k);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return toReturn;
	}

	double getLocalXeVolFracPoint(
			std::shared_ptr<xolotlSolver::PetscSolver> solver,
			int i, int j, int k) {
		double toReturn;
		try {
			// Get the solver handler
			auto& solverHandler = solver->getSolverHandler();
			// Get the rate at this location
			toReturn = solverHandler.getXeVolFracPoint(i, j, k);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return toReturn;
	}

	// std::vector<std::vector<std::vector<double> > > * getLocalXeConc(
	// 		std::shared_ptr<xolotlSolver::PetscSolver> solver) {
	// 	std::vector<std::vector<std::vector<double> > > * toReturn;
	// 	try {
	// 		// Get the solver handler
	// 		auto& solverHandler = solver->getSolverHandler();
	// 		// Get the rate at this location
	// 		toReturn = solverHandler.getLocalXeConc();
	// 	} catch (const std::exception& e) {
	// 		std::cerr << e.what() << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (const std::string& error) {
	// 		std::cerr << error << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (...) {
	// 		std::cerr << "Unrecognized exception seen." << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	}
	//
	// 	return toReturn;
	// }

	// double getLocalXeConc(
	// 		std::shared_ptr<xolotlSolver::PetscSolver> solver,
	// 		int i, int j, int k) {
	// 	double toReturn;
	// 	try {
	// 		// Get the solver handler
	// 		auto& solverHandler = solver->getSolverHandler();
	// 		// Get the rate at this location
	// 		toReturn = solverHandler.getLocalXeConc(i, j, k);
	// 	} catch (const std::exception& e) {
	// 		std::cerr << e.what() << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (const std::string& error) {
	// 		std::cerr << error << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (...) {
	// 		std::cerr << "Unrecognized exception seen." << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	}
	//
	// 	return toReturn;
	// }

	void getLocalCoordinates(
			std::shared_ptr<xolotlSolver::PetscSolver> solver, int *xs, int *xm,
			int *Mx, int *ys, int *ym, int *My, int *zs, int *zm, int *Mz) {
		try {
			// Get the solver handler
			auto& solverHandler = solver->getSolverHandler();
			// Get the local coordinates
			solverHandler.getLocalCoordinatesCpy(xs, xm, Mx, ys, ym, My, zs,
					zm, Mz);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return;
	}

	void initGBLocation(
			std::shared_ptr<xolotlSolver::PetscSolver> solver) {
		try {
			solver->initGBLocation();
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return;
	}

	void setGBLocation(
			std::shared_ptr<xolotlSolver::PetscSolver> solver, int i, int j,
			int k) {
		try {
			// Get the solver handler
			auto& solverHandler = solver->getSolverHandler();
			// Set the coordinate of the GB
			solverHandler.setGBLocation(i, j, k);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return;

	}

	// void setGBLocations(
	// 		std::shared_ptr<xolotlSolver::PetscSolver> solver,
	// 	  std::vector<std::tuple<int, int, int> > gbList) {
	// 	try {
	// 		// Get the solver handler
	// 		auto& solverHandler = solver->getSolverHandler();
	// 		// Set the coordinate of the GB
	// 		solverHandler.setGBLocations(gbList);
	// 	} catch (const std::exception& e) {
	// 		std::cerr << e.what() << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (const std::string& error) {
	// 		std::cerr << error << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (...) {
	// 		std::cerr << "Unrecognized exception seen." << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	}
	//
	// 	return;
	//
	// }

	void resetGBVector(std::shared_ptr<xolotlSolver::PetscSolver> solver) {
    try {
            // Get the solver handler
            auto& solverHandler = solver->getSolverHandler();
            // Reset the location
            solverHandler.resetGBVector();
    } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            std::cerr << "Aborting." << std::endl;
    } catch (const std::string& error) {
            std::cerr << error << std::endl;
            std::cerr << "Aborting." << std::endl;
    } catch (...) {
            std::cerr << "Unrecognized exception seen." << std::endl;
            std::cerr << "Aborting." << std::endl;
    }

    return;
	}

	void finalizeXolotl(
		std::shared_ptr<xolotlSolver::PetscSolver> solver, bool isStandalone) {
		try {
			// Call solver finalize
			//solver->finalize();
			solver->finalize(isStandalone);
		} catch (const std::exception& e) {
			std::cerr << e.what() << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (const std::string& error) {
			std::cerr << error << std::endl;
			std::cerr << "Aborting." << std::endl;
		} catch (...) {
			std::cerr << "Unrecognized exception seen." << std::endl;
			std::cerr << "Aborting." << std::endl;
		}

		return;
	}


	// TS getXolotlTS(
	// 	std::shared_ptr<xolotlSolver::PetscSolver> solver) {
	// 	TS toReturn = nullptr;
	// 	try {
	// 		// Set the time in the solver
	// 		solver->getTS(&toReturn);
	// 	} catch (const std::exception& e) {
	// 		std::cerr << e.what() << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (const std::string& error) {
	// 		std::cerr << error << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	} catch (...) {
	// 		std::cerr << "Unrecognized exception seen." << std::endl;
	// 		std::cerr << "Aborting." << std::endl;
	// 	}
	//
	// 	return toReturn;
	// }
};

extern "C" XolotlDLinterface* create(){
	return static_cast<XolotlDLinterface*>(new XolotlInterface);
}

extern "C" void destroy(XolotlDLinterface *ptr){
	delete ptr;
}
