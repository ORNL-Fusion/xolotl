#include <iostream>
#include <string>
#include "boost/program_options.hpp"
#include <XFile.h>
#include <MPIUtils.h>


namespace bpo = boost::program_options;
namespace xcore = xolotlCore;



int
main(int argc, char* argv[]) {

    int ret = 0;


    try {

        MPI_Init(&argc, &argv);


    	// Get the MPI communicator
    	auto xolotlComm = xolotlCore::MPIUtils::getMPIComm();
        // Determine our place in the MPI world.
        int cwRank;
        int cwSize;
        MPI_Comm_rank(xolotlComm, &cwRank);
        MPI_Comm_size(xolotlComm, &cwSize);

        // We have to be an MPI program,
        // but for now we only support a single process.
        if(cwSize > 1) {
            throw std::runtime_error("Multiple MPI processes not yet supported.  Run with 1 proc.");
        }

        // Parse the command line options.
        bool shouldRun = true;
        bpo::options_description desc("Supported options");
        desc.add_options()
            ("help", "show this help message")
            ("infile", bpo::value<std::string>(), "input file name")
        ;

        bpo::variables_map opts;
        bpo::store(bpo::parse_command_line(argc, argv, desc), opts);
        bpo::notify(opts);

        if(opts.count("help")) {
            std::cout << desc << '\n';
            shouldRun = false;
        }

        if((opts.count("infile") == 0) or
                opts["infile"].as<std::string>().empty()) {
            std::cerr << "input file name must not be empty" << std::endl;
            shouldRun = false;
            ret = 1;
        }

        if(shouldRun) {

            std::string fname = opts["infile"].as<std::string>();

            // Open the file.
            xcore::XFile xfile(fname,
            		xolotlComm,
                    xolotlCore::XFile::AccessMode::OpenReadWrite);

            // Determine the number of grid points.
            auto headerGroup = xfile.getGroup<xcore::XFile::HeaderGroup>();
            assert(headerGroup);
            xcore::HDF5File::Attribute<int> nxAttr(*headerGroup, "nx");
            auto nx = nxAttr.get();

            // Determine the last timestep written to the file.
            auto concGroup = xfile.getGroup<xcore::XFile::ConcentrationGroup>();
            assert(concGroup);
            xcore::HDF5File::Attribute<int> lastTimestepAttr(*concGroup, "lastTimeStep");
            auto lastTimeStep = lastTimestepAttr.get();


            std::cout << "nx: " << nx << '\n'
                << "last time step: " << lastTimeStep 
                << std::endl;

            // Open the timestep group associated with the 
            // last written timestep.
            auto tsGroup = concGroup->getLastTimestepGroup();
            assert(tsGroup);

            // Convert the last written timestep's concentrations to 
            // the new representation.
            xcore::XFile::TimestepGroup::Concs1DType allConcs(nx);
            for(auto x = 0; x < nx; ++x) {
                // Read the concentrations for the current position.
                std::ostringstream dsNameStr;
                dsNameStr << "position_" << x << "-1_-1";
                std::string dsName = dsNameStr.str();
                std::cout << "Reading conc data for gridpoint " << x << " from dataset " << dsName << std::endl;
                auto oldData = tsGroup->readGridPoint(x);
                auto nConcs = oldData.size();

                // Store into our ragged 2D representation.
                for(auto i = 0; i < nConcs; ++i) {
                    auto l = oldData[i][0];
                    auto conc = oldData[i][1];
                    allConcs[x].emplace_back(l, conc);
                }
            }

            // Write the dataset to the file.
            tsGroup->writeConcentrations(xfile, 0, allConcs);
        }
    }
    catch(std::exception& e) {
        std::cerr << e.what() << std::endl;
        ret = 1;
    }
    catch(...) {
        std::cerr << "Unrecognized exception caught." << std::endl;
        ret = 1;
    }

    // clean up
    MPI_Finalize();

    return ret;
}

