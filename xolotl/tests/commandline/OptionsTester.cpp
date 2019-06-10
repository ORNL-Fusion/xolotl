#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <Options.h>
#include <fstream>
#include <iostream>

using namespace std;

/**
 * Test suite for the Options class.
 */
BOOST_AUTO_TEST_SUITE (Options_testSuite)

BOOST_AUTO_TEST_CASE(noOptions) {
	xolotlCore::Options opts;

	// Build a fake, empty command line.
	// Here, "empty" is what the shell would give us if no arguments
	// were provided (i.e., it contains the executable name or path),
	// and we skipped that executable name before calling the
	// parsing method as specified in its comment.
	char* args[2];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = NULL;
	char** fargv = args;

	// Attempt to read the parameters
	fargv += 1;
	opts.readParams(fargv);

	// The Option class does not like empty command line
	// a parameter file is always needed
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);
}

BOOST_AUTO_TEST_CASE(badParamFileName) {
	xolotlCore::Options opts;

	string pathToFile("bla.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a non existing parameter file
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameters
	fargv += 1;
	opts.readParams(fargv);

	// If the parameter file does not exist, xolotl should not run
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);
}

BOOST_AUTO_TEST_CASE(badParamFile) {
	xolotlCore::Options opts;

	// Create a bad parameter file
	std::ofstream badParamFile("param_bad.txt");
	badParamFile << "netFile=tungsten.txt" << std::endl;
	badParamFile.close();

	string pathToFile("param_bad.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing bad options
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// An unrecognized parameter should result in indicating
	// the program shouldn't run, and an error exit code.
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);

	// Remove the created file
	std::string tempFile = "param_bad.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(goodParamFile) {
	xolotlCore::Options opts;

	// Create a good parameter file
	std::ofstream goodParamFile("param_good.txt");
	goodParamFile << "vizHandler=std" << std::endl
			<< "petscArgs=-fieldsplit_0_pc_type redundant "
					"-ts_max_snes_failures 200 "
					"-pc_fieldsplit_detect_coupling "
					"-ts_adapt_dt_max 10 -pc_type fieldsplit "
					"-fieldsplit_1_pc_type sor -ts_final_time 1000 "
					"-ts_max_steps 3" << std::endl << "networkFile=tungsten.txt"
			<< std::endl << "startTemp=900" << std::endl << "perfHandler=std"
			<< std::endl << "flux=1.5" << std::endl << "material=W100"
			<< std::endl << "initialV=0.05" << std::endl << "dimensions=1"
			<< std::endl << "voidPortion=60.0" << std::endl << "regularGrid=no"
			<< std::endl << "process=diff" << std::endl << "grouping=11 2 4"
			<< std::endl << "sputtering=0.5" << std::endl << "boundary=1 1"
			<< std::endl << "burstingDepth=5.0" << std::endl << "zeta=0.6"
			<< std::endl << "resoSize=10" << std::endl << "density=9.0"
			<< std::endl;
	goodParamFile.close();

	string pathToFile("param_good.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing good options
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should run with good parameters
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), true);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_SUCCESS);

	// Check the network filename
	BOOST_REQUIRE_EQUAL(opts.getNetworkFilename(), "tungsten.txt");

	// Check the temperature
	BOOST_REQUIRE_EQUAL(opts.useConstTemperatureHandlers(), true);
	BOOST_REQUIRE_EQUAL(opts.getConstTemperature(), 900.0);
	BOOST_REQUIRE_EQUAL(opts.getBulkTemperature(), 0.0);

	// Check if the flux option is used
	BOOST_REQUIRE_EQUAL(opts.useFluxAmplitude(), true);
	BOOST_REQUIRE_EQUAL(opts.getFluxAmplitude(), 1.5);

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerType(),
			xolotlPerf::IHandlerRegistry::std);

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.useVizStandardHandlers(), true);

	// Check the material option
	BOOST_REQUIRE_EQUAL(opts.getMaterial(), "W100");

	// Check the initial vacancy concentration option
	BOOST_REQUIRE_EQUAL(opts.getInitialVConcentration(), 0.05);

	// Check the number of dimensions option
	BOOST_REQUIRE_EQUAL(opts.getDimensionNumber(), 1);

	// Check the void portion option
	BOOST_REQUIRE_EQUAL(opts.getVoidPortion(), 60.0);

	// Check the regular grid option
	BOOST_REQUIRE_EQUAL(opts.useRegularXGrid(), false);

	// Check the grouping option
	BOOST_REQUIRE_EQUAL(opts.getGroupingMin(), 11);
	BOOST_REQUIRE_EQUAL(opts.getGroupingWidthA(), 2);
	BOOST_REQUIRE_EQUAL(opts.getGroupingWidthB(), 4);

	// Check the sputtering option
	BOOST_REQUIRE_EQUAL(opts.getSputteringYield(), 0.5);

	// Check the bursting depth option
	BOOST_REQUIRE_EQUAL(opts.getBurstingDepth(), 5.0);

	// Check the boundary conditions
	BOOST_REQUIRE_EQUAL(opts.getLeftBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getRightBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getBottomBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getTopBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getFrontBoundary(), 1);
	BOOST_REQUIRE_EQUAL(opts.getBackBoundary(), 1);

	// Check the electronic stopping power option
	BOOST_REQUIRE_EQUAL(opts.getZeta(), 0.6);

	// Check the minimum re-solution size option
	BOOST_REQUIRE_EQUAL(opts.getResoMinSize(), 10);

	// Check the xenon density option
	BOOST_REQUIRE_EQUAL(opts.getDensity(), 9.0);

	// Check the physical processes option
	auto map = opts.getProcesses();
	BOOST_REQUIRE_EQUAL(map["diff"], true);
	BOOST_REQUIRE_EQUAL(map["advec"], false);
	BOOST_REQUIRE_EQUAL(map["modifiedTM"], false);
	BOOST_REQUIRE_EQUAL(map["movingSurface"], false);
	BOOST_REQUIRE_EQUAL(map["reaction"], false);
	BOOST_REQUIRE_EQUAL(map["attenuation"], false);
	BOOST_REQUIRE_EQUAL(map["bursting"], false);

	// Check the PETSc options
	BOOST_REQUIRE_EQUAL(opts.getPetscArgc(), 16);

	// Get the detailed PETSc argument
	auto petscArgv = opts.getPetscArgv();
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[1], "-fieldsplit_0_pc_type"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[2], "redundant"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[3], "-ts_max_snes_failures"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[4], "200"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[5], "-pc_fieldsplit_detect_coupling"),
			0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[6], "-ts_adapt_dt_max"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[7], "10"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[8], "-pc_type"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[9], "fieldsplit"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[10], "-fieldsplit_1_pc_type"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[11], "sor"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[12], "-ts_final_time"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[13], "1000"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[14], "-ts_max_steps"), 0);
	BOOST_REQUIRE_EQUAL(strcmp(petscArgv[15], "3"), 0);

	// Remove the created file
	std::string tempFile = "param_good.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(goodParamFileNoHDF5) {
	xolotlCore::Options opts;

	// Create a good parameter file
	std::ofstream goodParamFile("param_good.txt");
	goodParamFile << "netParam=8 1 0 5 3" << std::endl << "grid=100 0.5"
			<< std::endl;
	goodParamFile.close();

	string pathToFile("param_good.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing good options
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should run with good parameters
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), true);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_SUCCESS);

	// Check if we use the HDF5 file
	BOOST_REQUIRE_EQUAL(opts.useHDF5(), false);

	// Check the network parameters
	BOOST_REQUIRE_EQUAL(opts.getMaxImpurity(), 8);
	BOOST_REQUIRE_EQUAL(opts.getMaxD(), 1);
	BOOST_REQUIRE_EQUAL(opts.getMaxT(), 0);
	BOOST_REQUIRE_EQUAL(opts.getMaxV(), 5);
	BOOST_REQUIRE_EQUAL(opts.getMaxI(), 3);
	BOOST_REQUIRE_EQUAL(opts.usePhaseCut(), false);

	// Check the grid parameters
	BOOST_REQUIRE_EQUAL(opts.getNX(), 100);
	BOOST_REQUIRE_EQUAL(opts.getXStepSize(), 0.5);
	BOOST_REQUIRE_EQUAL(opts.getNY(), 0);
	BOOST_REQUIRE_EQUAL(opts.getYStepSize(), 0.0);

	// Remove the created file
	std::string tempFile = "param_good.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(wrongPerfHandler) {
	xolotlCore::Options opts;

	// Create a parameter file with a wrong performance handler name
	std::ofstream paramFile("param_perf_wrong.txt");
	paramFile << "perfHandler=bogus" << std::endl;
	paramFile.close();

	string pathToFile("param_perf_wrong.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing a wrong performance handler option
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should not be able to run with a wrong performance handler parameter
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);

	// Remove the created file
	std::string tempFile = "param_perf_wrong.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(wrongVizHandler) {
	xolotlCore::Options opts;

	// Create a parameter file with a wrong visualization handler name
	std::ofstream paramFile("param_viz_wrong.txt");
	paramFile << "vizHandler=bogus" << std::endl;
	paramFile.close();

	string pathToFile("param_viz_wrong.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing a wrong performance handler option
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should not be able to run with a wrong visualization handler parameter
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);

	// Remove the created file
	std::string tempFile = "param_viz_wrong.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(goodParamFileWithProfiles) {
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

	xolotlCore::Options opts;

	// Create a parameter file using these two profile files
	std::ofstream paramFile("param_good_profiles.txt");
	paramFile << "fluxFile=fluxFile.dat" << std::endl
			<< "tempFile=temperatureFile.dat" << std::endl;
	paramFile.close();

	string pathToFile("param_good_profiles.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing
	// the two profile options
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should run with good parameters
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), true);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_SUCCESS);

	// Check the temperature
	BOOST_REQUIRE_EQUAL(opts.useTemperatureProfileHandlers(), true);
	BOOST_REQUIRE_EQUAL(opts.getTempProfileFilename(), "temperatureFile.dat");

	// Check if the heFlux option is used
	BOOST_REQUIRE_EQUAL(opts.useFluxAmplitude(), false);

	// Check if the time profile option is used for the flux
	BOOST_REQUIRE_EQUAL(opts.useFluxTimeProfile(), true);
	BOOST_REQUIRE_EQUAL(opts.getFluxProfileName(), "fluxFile.dat");

	// Remove the created files
	std::string tempFile = "temperatureFile.dat";
	std::remove(tempFile.c_str());
	tempFile = "fluxFile.dat";
	std::remove(tempFile.c_str());
	tempFile = "param_good_profiles.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(wrongFluxProfile) {
	xolotlCore::Options opts;

	// Create a parameter file with a wrong flux profile file name name
	std::ofstream paramFile("param_flux_wrong.txt");
	paramFile << "fluxFile=bogus" << std::endl;
	paramFile.close();

	string pathToFile("param_flux_wrong.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing a wrong performance handler option
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should not be able to run with a wrong flux profile file name
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);

	// Remove the created file
	std::string tempFile = "param_flux_wrong.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(wrongTempProfile) {
	xolotlCore::Options opts;

	// Create a parameter file with a wrong temperature profile file name name
	std::ofstream paramFile("param_temp_wrong.txt");
	paramFile << "tempFile=bogus" << std::endl;
	paramFile.close();

	string pathToFile("param_temp_wrong.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing a wrong performance handler option
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should not be able to run with a wrong temperature profile file name
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), false);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_FAILURE);

	// Remove the created file
	std::string tempFile = "param_temp_wrong.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(papiPerfHandler) {
	xolotlCore::Options opts;

	// Create a parameter file using the PAPI performance handlers
	std::ofstream paramFile("param_good_perf_papi.txt");
	paramFile << "perfHandler=papi" << std::endl;
	paramFile.close();

	string pathToFile("param_good_perf_papi.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should run with good parameters
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), true);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_SUCCESS);

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerType(),
			xolotlPerf::IHandlerRegistry::papi);

	// Remove the created file
	std::string tempFile = "param_good_perf_papi.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(osPerfHandler) {
	xolotlCore::Options opts;

	// Create a parameter file using the OS performance handlers
	std::ofstream paramFile("param_good_perf_os.txt");
	paramFile << "perfHandler=os" << std::endl;
	paramFile.close();

	string pathToFile("param_good_perf_os.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should run with good parameters
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), true);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_SUCCESS);

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerType(),
			xolotlPerf::IHandlerRegistry::os);

	// Remove the created file
	std::string tempFile = "param_good_perf_os.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_CASE(dummyPerfHandler) {
	xolotlCore::Options opts;

	// Create a parameter file using the dummy performance handlers
	std::ofstream paramFile("param_good_perf_dummy.txt");
	paramFile << "perfHandler=dummy" << std::endl;
	paramFile.close();

	string pathToFile("param_good_perf_dummy.txt");
	string filename = pathToFile;
	const char* fname = filename.c_str();

	// Build a command line with a parameter file containing good options
	char* args[3];
	args[0] = const_cast<char*>("./xolotl");
	args[1] = const_cast<char*>(fname);
	args[2] = NULL;
	char** fargv = args;

	// Attempt to read the parameter file
	fargv += 1;
	opts.readParams(fargv);

	// Xolotl should run with good parameters
	BOOST_REQUIRE_EQUAL(opts.shouldRun(), true);
	BOOST_REQUIRE_EQUAL(opts.getExitCode(), EXIT_SUCCESS);

	// Check the performance handler
	BOOST_REQUIRE_EQUAL(opts.getPerfHandlerType(),
			xolotlPerf::IHandlerRegistry::dummy);

	// Remove the created file
	std::string tempFile = "param_good_perf_dummy.txt";
	std::remove(tempFile.c_str());
}

BOOST_AUTO_TEST_SUITE_END()
