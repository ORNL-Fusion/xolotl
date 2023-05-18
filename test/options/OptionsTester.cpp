#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <fstream>
#include <iostream>

#include <boost/program_options/errors.hpp>
#include <boost/test/framework.hpp>
#include <boost/test/unit_test.hpp>
namespace bpo = boost::program_options;

#include <xolotl/options/Options.h>

using namespace std;
using namespace xolotl::options;

/**
 * Test suite for the Options class.
 */
BOOST_AUTO_TEST_SUITE(Options_testSuite)

BOOST_AUTO_TEST_CASE(noOptions)
{
	Options opts;

	// Build a fake, empty command line.
	// Here, "empty" is what the shell would give us if no arguments
	// were provided (i.e., it contains the executable name or path),
	// and we skipped that executable name before calling the
	// parsing method as specified in its comment.
	const char* argv[] = {"./xolotl"};

	// Attempt to read the parameters
	BOOST_CHECK_THROW(opts.readParams(1, argv), bpo::error);
}

BOOST_AUTO_TEST_CASE(badParamFileName)
{
	Options opts;

	string pathToFile("bla.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a non existing parameter file
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameters
	BOOST_CHECK_THROW(opts.readParams(2, argv), bpo::error);
}

BOOST_AUTO_TEST_CASE(badParamFile)
{
	try {
		Options opts;

		// Create a bad parameter file
		std::ofstream badParamFile("param_bad.txt");
		badParamFile << "netFile=tungsten.txt" << std::endl;
		badParamFile.close();

		string pathToFile("param_bad.txt");
		string filename = pathToFile;
		const char* fname = filename.c_str();

		// Build a command line with a parameter file containing bad options
		const char* argv[] = {"./xolotl", fname};

		// Attempt to read the parameter file
		opts.readParams(2, argv);

		// Remove the created file
		std::string tempFile = "param_bad.txt";
		std::remove(tempFile.c_str());
	}
	catch (const std::exception& e) {
		// Great
		std::cerr << e.what() << std::endl;
		return;
	}

	// Should not get here
	BOOST_REQUIRE(false);
}

BOOST_AUTO_TEST_CASE(goodParamFile)
{
	Options opts;

	// Create a good parameter file
	std::ofstream goodParamFile("param_good.txt");
	goodParamFile
		<< "vizHandler=std" << std::endl
		<< "petscArgs=-fieldsplit_0_pc_type redundant "
		   "-ts_max_snes_failures 200 "
		   "-pc_fieldsplit_detect_coupling "
		   "-ts_adapt_dt_max 10 -pc_type fieldsplit "
		   "-fieldsplit_1_pc_type sor -ts_final_time 1000 "
		   "-ts_max_steps 3"
		<< std::endl
		<< "networkFile=tungsten.txt" << std::endl
		<< "tempHandler=constant" << std::endl
		<< "tempParam=900" << std::endl
		<< "perfHandler=os" << std::endl
		<< "flux=1.5" << std::endl
		<< "material=W100" << std::endl
		<< "initialConc=V 1 0.05" << std::endl
		<< "dimensions=1" << std::endl
		<< "gridType=nonuniform" << std::endl
		<< "gridParam=10" << std::endl
		<< "process=diff" << std::endl
		<< "grouping=11 2 4" << std::endl
		<< "sputtering=0.5" << std::endl
		<< "boundary=1 1" << std::endl
		<< "burstingDepth=5.0" << std::endl
		<< "burstingFactor=2.5" << std::endl
		<< "zeta=0.6" << std::endl
		<< "radiusSize=5 0 3" << std::endl
		<< "density=9.0" << std::endl
		<< "lattice=0.1" << std::endl
		<< "impurityRadius=0.5" << std::endl
		<< "biasFactor=2.0" << std::endl
		<< "hydrogenFactor=0.5" << std::endl
		<< "xenonDiffusivity=3.0" << std::endl
		<< "fissionYield=0.3" << std::endl
		<< "heVRatio=5.0" << std::endl
		<< "migrationThreshold=1.0" << std::endl
		<< "fluxDepthProfileFilePath=path/to/the/flux/profile/file.txt"
		<< std::endl;
	goodParamFile.close();

	string pathToFile("param_good.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing good options
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the network filename
	BOOST_REQUIRE_EQUAL(opts.getNetworkFilename(), "tungsten.txt");

	// Check the temperature
	BOOST_REQUIRE_EQUAL(opts.getTempHandlerName(), "constant");
	BOOST_REQUIRE_EQUAL(opts.getTempParam(), 900.0);
	BOOST_REQUIRE_EQUAL(opts.getTempParam(1), 0.0);

	// Check if the flux option is used
	BOOST_REQUIRE_EQUAL(opts.useFluxAmplitude(), true);
	BOOST_REQUIRE_EQUAL(opts.getFluxAmplitude(), 1.5);

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerName(), "os");

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getVizHandlerName(), "std");

	// Check the material option
	BOOST_REQUIRE_EQUAL(opts.getMaterial(), "W100");

	// Check the initial vacancy concentration option
	BOOST_REQUIRE_EQUAL(opts.getInitialConcentration(), "V 1 0.05");

	// Check the number of dimensions option
	BOOST_REQUIRE_EQUAL(opts.getDimensionNumber(), 1);

	// Check the grid options
	BOOST_REQUIRE_EQUAL(opts.getGridTypeName(), "nonuniform");
	BOOST_REQUIRE_EQUAL(opts.getGridParam(0), 10);

	// Check the grouping option
	BOOST_REQUIRE_EQUAL(opts.getGroupingMin(), 11);
	BOOST_REQUIRE_EQUAL(opts.getGroupingWidthA(), 2);
	BOOST_REQUIRE_EQUAL(opts.getGroupingWidthB(), 4);

	// Check the sputtering option
	BOOST_REQUIRE_EQUAL(opts.getSputteringYield(), 0.5);

	// Check the bursting depth option
	BOOST_REQUIRE_EQUAL(opts.getBurstingDepth(), 5.0);

	// Check the bursting factor option
	BOOST_REQUIRE_EQUAL(opts.getBurstingFactor(), 2.5);

	// Check the boundary conditions
	BOOST_REQUIRE_EQUAL(opts.getLeftBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getRightBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getBottomBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getTopBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getFrontBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getBackBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getBCString(), "mirror");

	// Check the electronic stopping power option
	BOOST_REQUIRE_EQUAL(opts.getZeta(), 0.6);

	// Check the minimum size option for the average radius computation
	auto sizes = opts.getRadiusMinSizes();
	BOOST_REQUIRE_EQUAL(sizes[0], 5);
	BOOST_REQUIRE_EQUAL(sizes[1], 0);
	BOOST_REQUIRE_EQUAL(sizes[2], 3);

	// Check the xenon density option
	BOOST_REQUIRE_EQUAL(opts.getDensity(), 9.0);

	// Check the lattice parameter option
	BOOST_REQUIRE_EQUAL(opts.getLatticeParameter(), 0.1);

	// Check the impurity radius option
	BOOST_REQUIRE_EQUAL(opts.getImpurityRadius(), 0.5);

	// Check the interstitial bias factor option
	BOOST_REQUIRE_EQUAL(opts.getBiasFactor(), 2.0);

	// Check the hydrogen radius factor option
	BOOST_REQUIRE_EQUAL(opts.getHydrogenFactor(), 0.5);

	// Check the xenon diffusivity option
	BOOST_REQUIRE_EQUAL(opts.getXenonDiffusivity(), 3.0);

	// Check the fission yield option
	BOOST_REQUIRE_EQUAL(opts.getFissionYield(), 0.3);

	// Check the HeV ration option
	BOOST_REQUIRE_EQUAL(opts.getHeVRatio(), 5.0);

	// Check the migration threshold option
	BOOST_REQUIRE_EQUAL(opts.getMigrationThreshold(), 1.0);

	// Check the network filename
	BOOST_REQUIRE_EQUAL(opts.getFluxDepthProfileFilePath(),
		"path/to/the/flux/profile/file.txt");

	// Check the physical processes option
	auto map = opts.getProcesses();
	BOOST_REQUIRE_EQUAL(map["diff"], true);
	BOOST_REQUIRE_EQUAL(map["advec"], false);
	BOOST_REQUIRE_EQUAL(map["modifiedTM"], false);
	BOOST_REQUIRE_EQUAL(map["movingSurface"], false);
	BOOST_REQUIRE_EQUAL(map["reaction"], false);
	BOOST_REQUIRE_EQUAL(map["attenuation"], false);
	BOOST_REQUIRE_EQUAL(map["bursting"], false);
	BOOST_REQUIRE_EQUAL(map["resolution"], false);
	BOOST_REQUIRE_EQUAL(map["heterogeneous"], false);

	// Check the PETSc options
	BOOST_REQUIRE_EQUAL(opts.getPetscArg(),
		"-fieldsplit_0_pc_type redundant "
		"-ts_max_snes_failures 200 "
		"-pc_fieldsplit_detect_coupling "
		"-ts_adapt_dt_max 10 -pc_type fieldsplit "
		"-fieldsplit_1_pc_type sor -ts_final_time 1000 "
		"-ts_max_steps 3");

	// Remove the created file
	std::string tempFile = "param_good.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(goodParamFileNoHDF5)
{
	Options opts;

	// Create a good parameter file
	std::ofstream goodParamFile("param_good.txt");
	goodParamFile << "netParam=8 1 0 5 3" << std::endl;
	goodParamFile.close();

	string pathToFile("param_good.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing good options
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check if we use the HDF5 file
	BOOST_REQUIRE_EQUAL(opts.useHDF5(), false);

	// Check the network parameters
	BOOST_REQUIRE_EQUAL(opts.getMaxImpurity(), 8);
	BOOST_REQUIRE_EQUAL(opts.getMaxD(), 1);
	BOOST_REQUIRE_EQUAL(opts.getMaxT(), 0);
	BOOST_REQUIRE_EQUAL(opts.getMaxV(), 5);
	BOOST_REQUIRE_EQUAL(opts.getMaxI(), 3);

	// Remove the created file
	std::string tempFile = "param_good.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(wrongPerfHandler)
{
	Options opts;

	// Create a parameter file with a wrong performance handler name
	std::ofstream paramFile("param_perf_wrong.txt");
	paramFile << "perfHandler=bogus" << std::endl;
	paramFile.close();

	string pathToFile("param_perf_wrong.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_CHECK_THROW(opts.readParams(2, argv), bpo::invalid_option_value);

	// Remove the created file
	std::string tempFile = "param_perf_wrong.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(wrongVizHandler)
{
	Options opts;

	// Create a parameter file with a wrong visualization handler name
	std::ofstream paramFile("param_viz_wrong.txt");
	paramFile << "vizHandler=bogus" << std::endl;
	paramFile.close();

	string pathToFile("param_viz_wrong.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_CHECK_THROW(opts.readParams(2, argv), bpo::invalid_option_value);

	// Remove the created file
	std::string tempFile = "param_viz_wrong.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(goodParamFileWithProfiles)
{
	// Create a file with temperature profile data
	// First column with the time and the second with
	// the temperature at that time.
	std::ofstream writeTempFile("temperatureFile.dat");
	writeTempFile << "0.0 2.0 \n"
					 "1.0 1.99219766723 \n"
					 "2.0 1.87758256189 \n"
					 "3.0 1.4311765168 \n"
					 "4.0 0.583853163453 \n"
					 "5.0 0.000137654918313 \n"
					 "6.0 0.789204200569 \n"
					 "7.0 1.9875147713 \n"
					 "8.0 0.854499966191 \n"
					 "9.0 0.235300873168 \n"
					 "10.0 1.99779827918";
	writeTempFile.close();

	// Create a file with a time profile for the flux
	// First column with the time and the second with
	// the amplitude (in He/nm2/s) at that time.
	std::ofstream writeFluxFile("fluxFile.dat");
	writeFluxFile << "0.0 1000.0 \n"
					 "1.0 4000.0 \n"
					 "2.0 2000.0 \n"
					 "3.0 3000.0 \n"
					 "4.0 0.0";
	writeFluxFile.close();

	Options opts;

	// Create a parameter file using these two profile files
	std::ofstream paramFile("param_good_profiles.txt");
	paramFile << "fluxFile=fluxFile.dat" << std::endl
			  << "tempHandler=profile" << std::endl
			  << "tempFile=temperatureFile.dat" << std::endl;
	paramFile.close();

	string pathToFile("param_good_profiles.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing
	// the two profile options
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the temperature
	BOOST_REQUIRE_EQUAL(opts.getTempHandlerName(), "profile");
	BOOST_REQUIRE_EQUAL(opts.getTempProfileFilename(), "temperatureFile.dat");

	// Check if the heFlux option is used
	BOOST_REQUIRE_EQUAL(opts.useFluxAmplitude(), false);

	// Check if the time profile option is used for the flux
	BOOST_REQUIRE_EQUAL(opts.useFluxTimeProfile(), true);
	BOOST_REQUIRE_EQUAL(opts.getFluxTimeProfileFilePath(), "fluxFile.dat");

	// Remove the created files
	std::string tempFile = "temperatureFile.dat";
	std::remove(tempFile.c_str());
	tempFile = "fluxFile.dat";
	std::remove(tempFile.c_str());
	tempFile = "param_good_profiles.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(readGridIn)
{
	// Create a grid file
	std::ofstream gridFile("grid.dat");
	gridFile << "0.0 1.0 3.0 7.0";
	gridFile.close();

	Options opts;

	// Create a parameter file with a file name for the grid
	std::ofstream paramFile("param_read_in_grid.txt");
	paramFile << "gridType=read" << std::endl
			  << "gridFile=grid.dat" << std::endl;
	paramFile.close();

	string pathToFile("param_read_in_grid.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the grid options
	BOOST_REQUIRE_EQUAL(opts.getGridTypeName(), "read");
	BOOST_REQUIRE_EQUAL(opts.getGridFilename(), "grid.dat");

	// Remove the created file
	std::string tempFile = "param_read_in_grid.txt";
	std::remove(tempFile.c_str());
	tempFile = "grid.dat";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(wrongFluxProfile)
{
	Options opts;

	// Create a parameter file with a wrong flux profile file name name
	std::ofstream paramFile("param_flux_wrong.txt");
	paramFile << "fluxFile=bogus" << std::endl;
	paramFile.close();

	string pathToFile("param_flux_wrong.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_CHECK_THROW(opts.readParams(2, argv), bpo::error);

	// Remove the created file
	std::string tempFile = "param_flux_wrong.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(wrongTempProfile)
{
	Options opts;

	// Create a parameter file with a wrong temperature profile file name name
	std::ofstream paramFile("param_temp_wrong.txt");
	paramFile << "tempFile=bogus" << std::endl;
	paramFile.close();

	string pathToFile("param_temp_wrong.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing a wrong performance
	// handler option
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_CHECK_THROW(opts.readParams(2, argv), bpo::error);

	// Remove the created file
	std::string tempFile = "param_temp_wrong.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(papiPerfHandler)
{
	Options opts;

	// Create a parameter file using the PAPI performance handlers
	std::ofstream paramFile("param_good_perf_papi.txt");
	paramFile << "perfHandler=papi" << std::endl;
	paramFile.close();

	string pathToFile("param_good_perf_papi.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerName(), "papi");

	// Remove the created file
	std::string tempFile = "param_good_perf_papi.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(osPerfHandler)
{
	Options opts;

	// Create a parameter file using the OS performance handlers
	std::ofstream paramFile("param_good_perf_os.txt");
	paramFile << "perfHandler=os" << std::endl;
	paramFile.close();

	string pathToFile("param_good_perf_os.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerName(), "os");

	// Remove the created file
	std::string tempFile = "param_good_perf_os.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(dummyPerfHandler)
{
	Options opts;

	// Create a parameter file using the dummy performance handlers
	std::ofstream paramFile("param_good_perf_dummy.txt");
	paramFile << "perfHandler=dummy" << std::endl;
	paramFile.close();

	string pathToFile("param_good_perf_dummy.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing good options
	const char* argv[] = {"./xolotl", fname};

	// Attempt to read the parameter file
	BOOST_REQUIRE_NO_THROW(opts.readParams(2, argv));

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerName(), "dummy");

	// Remove the created file
	std::string tempFile = "param_good_perf_dummy.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
