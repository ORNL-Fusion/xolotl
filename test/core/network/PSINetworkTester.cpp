#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
namespace utf = boost::unit_test;

#include <xolotl/core/network/PSINetworkHandler.h>
#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/options/ConfOptions.h>
#include <xolotl/test/CommandLine.h>
#include <xolotl/test/MPITestUtils.h>
#include <xolotl/test/PSINetworkTesterData.h>

using namespace std;
using namespace xolotl;
using namespace core;
using namespace network;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

// Initialize MPI before running any tests; finalize it running all tests.
BOOST_GLOBAL_FIXTURE(MPIFixture);

/**
 * This suite is responsible for testing the PSI network.
 */
BOOST_AUTO_TEST_SUITE(PSINetwork)

BOOST_AUTO_TEST_CASE(fullyRefined)
{
	loadNetworkHandlers();
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 1 1 1 1" << std::endl
			  << "process=reaction" << std::endl
			  << "tempParam=1000" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIFullSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the temperature and lattice constant from the options
	double latticeConst = opts.getLatticeParameter() <= 0.0 ?
		xolotl::core::tungstenLatticeConstant :
		opts.getLatticeParameter();
	double temperature = opts.getTempParam();

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe =
		util::getMaxHePerVLoop(maxV, latticeConst, temperature);
	NetworkType::AmountType maxD = 2.0 / 3.0 * (double)maxHe;
	NetworkType::AmountType maxT = 2.0 / 3.0 * (double)maxHe;
	NetworkType network({maxHe, maxD, maxT, maxV, maxI}, 1, opts);

	BOOST_REQUIRE(network.hasDeuterium());
	BOOST_REQUIRE(network.hasTritium());

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 156);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 156);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.317, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.0159275, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 156);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 5);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 5);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 4);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 1, 3, 10, 32, 36, 40, 47, 58, 69, 85, 107, 31, 106, 35,
		84, 39, 68, 46, 57};
	knownDFill[1] = {1, 0, 2, 9, 31, 35, 39, 46, 57, 68, 84, 106, 3, 10, 32, 36,
		40, 47, 58, 69, 85, 107};
	knownDFill[2] = {2, 0, 3, 1, 4, 5, 6, 7, 10, 11, 12, 13, 14, 16, 17, 18, 19,
		21, 22, 23, 25, 26, 28, 32, 36, 40, 41, 43, 47, 48, 49, 51, 52, 54, 58,
		59, 60, 62, 63, 65, 69, 70, 71, 72, 74, 75, 76, 78, 79, 81, 85, 86, 87,
		88, 89, 91, 92, 93, 94, 96, 97, 98, 100, 101, 103, 107, 108, 109, 110,
		111, 113, 114, 115, 116, 118, 119, 120, 122, 123, 125, 128, 129, 130,
		131, 132, 133, 135, 136, 137, 138, 139, 141, 142, 143, 144, 146, 147,
		148, 150, 151, 153, 8, 15, 20, 24, 27, 29, 33, 37, 42, 44, 50, 53, 55,
		61, 64, 66, 73, 77, 80, 82, 90, 95, 99, 102, 104, 112, 117, 121, 124,
		126, 134, 140, 145, 149, 152, 154};
	knownDFill[3] = {3, 0, 1, 2, 9, 31, 35, 39, 46, 57, 68, 84, 4, 11, 33};
	knownDFill[4] = {4, 2, 3, 9, 39, 46, 57, 68, 84, 5, 12};
	knownDFill[5] = {5, 2, 4, 9, 46, 57, 68, 84, 6, 13};
	knownDFill[6] = {6, 2, 5, 9, 68, 84, 7, 14};
	knownDFill[7] = {7, 2, 6, 9, 84, 8, 15};
	knownDFill[8] = {8, 2, 7};
	knownDFill[9] = {9, 0, 10, 1, 3, 4, 5, 6, 7, 11, 12, 13, 14, 16, 17, 18, 19,
		21, 22, 23, 25, 26, 28, 32, 36, 40, 41, 43, 47, 48, 49, 51, 52, 54, 58,
		59, 60, 62, 63, 65, 69, 70, 71, 72, 74, 75, 76, 78, 79, 81, 85, 86, 87,
		88, 89, 91, 92, 93, 94, 96, 97, 98, 100, 101, 103, 107, 108, 109, 110,
		111, 113, 114, 115, 116, 118, 119, 120, 122, 123, 125, 128, 129, 130,
		131, 132, 133, 135, 136, 137, 138, 139, 141, 142, 143, 144, 146, 147,
		148, 150, 151, 153, 15, 20, 24, 27, 29, 30, 34, 38, 44, 45, 53, 55, 56,
		64, 66, 67, 77, 80, 82, 83, 95, 99, 102, 104, 105, 117, 121, 124, 126,
		127, 140, 145, 149, 152, 154, 155};
	knownDFill[10] = {10, 0, 1, 9, 2, 31, 35, 39, 46, 57, 68, 84, 11, 16, 34};
	knownDFill[11] = {11, 2, 10, 3, 9, 39, 46, 57, 68, 84, 12, 17};
	knownDFill[12] = {12, 2, 11, 4, 9, 46, 57, 68, 84, 13, 18};
	knownDFill[13] = {13, 2, 12, 5, 9, 68, 84, 14, 19};
	knownDFill[14] = {14, 2, 13, 6, 9, 84, 15, 20};
	knownDFill[15] = {15, 2, 14, 7, 9};
	knownDFill[16] = {16, 2, 9, 10, 39, 46, 57, 68, 84, 17, 21};
	knownDFill[17] = {17, 2, 16, 9, 11, 46, 57, 68, 84, 18, 22};
	knownDFill[18] = {18, 2, 17, 9, 12, 68, 84, 19, 23};
	knownDFill[19] = {19, 2, 18, 9, 13, 84, 20, 24};
	knownDFill[20] = {20, 2, 19, 9, 14};
	knownDFill[21] = {21, 2, 9, 16, 46, 57, 68, 84, 22, 25};
	knownDFill[22] = {22, 2, 21, 9, 17, 68, 84, 23, 26};
	knownDFill[23] = {23, 2, 22, 9, 18, 84, 24, 27};
	knownDFill[24] = {24, 2, 23, 9, 19};
	knownDFill[25] = {25, 2, 9, 21, 68, 84, 26, 28};
	knownDFill[26] = {26, 2, 25, 9, 22, 84, 27, 29};
	knownDFill[27] = {27, 2, 26, 9, 23};
	knownDFill[28] = {28, 2, 9, 25, 84, 29, 30};
	knownDFill[29] = {29, 2, 28, 9, 26};
	knownDFill[30] = {30, 9, 28};
	knownDFill[31] = {31, 0, 32, 1, 3, 10, 33, 34, 35, 36, 37, 38, 39, 40, 41,
		42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
		60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
		78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
		96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
		111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124,
		125, 126, 127, 128, 129, 130, 131, 132, 133, 135, 136, 137, 138, 139,
		141, 142, 143, 144, 146, 147, 148, 150, 151, 153};
	knownDFill[32] = {32, 0, 1, 31, 2, 9, 35, 39, 46, 57, 68, 84, 33, 34, 36};
	knownDFill[33] = {33, 2, 32, 3, 31, 35, 39, 46, 57, 68, 84, 37};
	knownDFill[34] = {34, 9, 32, 10, 31, 35, 39, 46, 57, 68, 84, 38};
	knownDFill[35] = {35, 0, 36, 1, 3, 10, 31, 32, 33, 34, 37, 38, 39, 40, 41,
		42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59,
		60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,
		78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
		96, 97, 98, 99, 100, 101, 102, 103, 104, 105};
	knownDFill[36] = {
		36, 0, 1, 35, 2, 9, 31, 32, 39, 46, 57, 68, 84, 37, 38, 40};
	knownDFill[37] = {37, 2, 36, 3, 35, 31, 33, 39, 46, 57, 68, 84, 41};
	knownDFill[38] = {38, 9, 36, 10, 35, 31, 34, 39, 46, 57, 68, 84, 43};
	knownDFill[39] = {39, 0, 40, 1, 3, 4, 10, 11, 16, 31, 35, 32, 33, 34, 36,
		37, 38, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56,
		57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74,
		75, 76, 77, 78, 79, 80, 81, 82, 83};
	knownDFill[40] = {
		40, 0, 1, 39, 2, 9, 31, 36, 32, 35, 46, 57, 68, 41, 43, 47};
	knownDFill[41] = {
		41, 2, 40, 3, 39, 9, 31, 37, 33, 35, 46, 57, 68, 42, 44, 48};
	knownDFill[42] = {42, 2, 41, 4, 39, 31, 35, 46, 57, 68, 49};
	knownDFill[43] = {
		43, 2, 9, 40, 10, 39, 31, 38, 34, 35, 46, 57, 68, 44, 45, 51};
	knownDFill[44] = {44, 2, 43, 9, 41, 11, 39, 31, 35, 46, 57, 68, 52};
	knownDFill[45] = {45, 9, 43, 16, 39, 31, 35, 46, 57, 68, 54};
	knownDFill[46] = {46, 0, 47, 1, 3, 4, 5, 10, 11, 12, 16, 17, 21, 31, 39, 32,
		33, 34, 35, 36, 37, 38, 40, 41, 42, 43, 44, 45, 48, 49, 50, 51, 52, 53,
		54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67};
	knownDFill[47] = {
		47, 0, 1, 46, 2, 9, 31, 40, 32, 39, 35, 36, 57, 48, 51, 58};
	knownDFill[48] = {
		48, 2, 47, 3, 46, 9, 31, 41, 33, 39, 35, 37, 57, 49, 52, 59};
	knownDFill[49] = {49, 2, 48, 4, 46, 9, 31, 42, 35, 39, 57, 50, 53, 60};
	knownDFill[50] = {50, 2, 49, 5, 46, 31, 35, 39, 57, 61};
	knownDFill[51] = {
		51, 2, 9, 47, 10, 46, 31, 43, 34, 39, 35, 38, 57, 52, 54, 62};
	knownDFill[52] = {52, 2, 51, 9, 48, 11, 46, 31, 44, 35, 39, 57, 53, 55, 63};
	knownDFill[53] = {53, 2, 52, 9, 49, 12, 46, 31, 35, 39, 57, 64};
	knownDFill[54] = {54, 2, 9, 51, 16, 46, 31, 45, 35, 39, 57, 55, 56, 65};
	knownDFill[55] = {55, 2, 54, 9, 52, 17, 46, 31, 35, 39, 57, 66};
	knownDFill[56] = {56, 9, 54, 21, 46, 31, 35, 39, 57, 67};
	knownDFill[57] = {57, 0, 58, 1, 3, 4, 5, 10, 11, 12, 16, 17, 21, 31, 46, 32,
		33, 34, 35, 39, 36, 37, 38, 40, 41, 42, 43, 44, 45, 47, 48, 49, 50, 51,
		52, 53, 54, 55, 56, 68};
	knownDFill[58] = {
		58, 0, 1, 57, 2, 9, 31, 47, 32, 46, 35, 40, 36, 39, 59, 62, 69};
	knownDFill[59] = {
		59, 2, 58, 3, 57, 9, 31, 48, 33, 46, 35, 41, 37, 39, 60, 63, 70};
	knownDFill[60] = {60, 2, 59, 4, 57, 9, 31, 49, 35, 42, 39, 46, 61, 64, 71};
	knownDFill[61] = {61, 2, 60, 5, 57, 31, 50, 35, 39, 46, 72};
	knownDFill[62] = {
		62, 2, 9, 58, 10, 57, 31, 51, 34, 46, 35, 43, 38, 39, 63, 65, 74};
	knownDFill[63] = {
		63, 2, 62, 9, 59, 11, 57, 31, 52, 35, 44, 39, 46, 64, 66, 75};
	knownDFill[64] = {64, 2, 63, 9, 60, 12, 57, 31, 53, 35, 39, 46, 76};
	knownDFill[65] = {65, 2, 9, 62, 16, 57, 31, 54, 35, 45, 39, 46, 66, 67, 78};
	knownDFill[66] = {66, 2, 65, 9, 63, 17, 57, 31, 55, 35, 39, 46, 79};
	knownDFill[67] = {67, 9, 65, 21, 57, 31, 56, 35, 39, 46, 81};
	knownDFill[68] = {68, 0, 69, 1, 3, 4, 5, 6, 10, 11, 12, 13, 16, 17, 18, 21,
		22, 25, 31, 57, 32, 33, 34, 35, 46, 36, 37, 38, 39, 40, 41, 42, 43, 44,
		45, 84};
	knownDFill[69] = {
		69, 0, 1, 68, 2, 9, 31, 58, 32, 57, 35, 47, 36, 46, 39, 40, 70, 74, 85};
	knownDFill[70] = {70, 2, 69, 3, 68, 9, 31, 59, 33, 57, 35, 48, 37, 46, 39,
		41, 71, 75, 86};
	knownDFill[71] = {71, 2, 70, 4, 68, 9, 31, 60, 35, 49, 39, 42, 72, 76, 87};
	knownDFill[72] = {72, 2, 71, 5, 68, 9, 31, 61, 35, 50, 39, 73, 77, 88};
	knownDFill[73] = {73, 2, 72, 6, 68, 31, 35, 39, 89};
	knownDFill[74] = {74, 2, 9, 69, 10, 68, 31, 62, 34, 57, 35, 51, 38, 46, 39,
		43, 75, 78, 91};
	knownDFill[75] = {
		75, 2, 74, 9, 70, 11, 68, 31, 63, 35, 52, 39, 44, 76, 79, 92};
	knownDFill[76] = {76, 2, 75, 9, 71, 12, 68, 31, 64, 35, 53, 39, 77, 80, 93};
	knownDFill[77] = {77, 2, 76, 9, 72, 13, 68, 31, 35, 39, 94};
	knownDFill[78] = {78, 2, 9, 74, 16, 68, 31, 65, 35, 54, 39, 45, 79, 81, 96};
	knownDFill[79] = {79, 2, 78, 9, 75, 17, 68, 31, 66, 35, 55, 39, 80, 82, 97};
	knownDFill[80] = {80, 2, 79, 9, 76, 18, 68, 31, 35, 39, 98};
	knownDFill[81] = {81, 2, 9, 78, 21, 68, 31, 67, 35, 56, 39, 82, 83, 100};
	knownDFill[82] = {82, 2, 81, 9, 79, 22, 68, 31, 35, 39, 101};
	knownDFill[83] = {83, 9, 81, 25, 68, 31, 35, 39, 103};
	knownDFill[84] = {84, 0, 85, 1, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 16, 17,
		18, 19, 21, 22, 23, 25, 26, 28, 31, 68, 32, 33, 34, 35, 57, 36, 37, 38,
		39, 46, 106};
	knownDFill[85] = {85, 0, 1, 84, 2, 9, 31, 69, 32, 68, 35, 58, 36, 57, 39,
		47, 40, 46, 86, 91, 107};
	knownDFill[86] = {86, 2, 85, 3, 84, 9, 31, 70, 33, 68, 35, 59, 37, 57, 39,
		48, 41, 46, 87, 92, 108};
	knownDFill[87] = {
		87, 2, 86, 4, 84, 9, 31, 71, 35, 60, 39, 49, 42, 46, 88, 93, 109};
	knownDFill[88] = {88, 2, 87, 5, 84, 9, 31, 72, 35, 61, 39, 50, 89, 94, 110};
	knownDFill[89] = {89, 2, 88, 6, 84, 9, 31, 73, 35, 90, 95, 111};
	knownDFill[90] = {90, 2, 89, 7, 84, 31, 35, 112};
	knownDFill[91] = {91, 2, 9, 85, 10, 84, 31, 74, 34, 68, 35, 62, 38, 57, 39,
		51, 43, 46, 92, 96, 113};
	knownDFill[92] = {
		92, 2, 91, 9, 86, 11, 84, 31, 75, 35, 63, 39, 52, 44, 46, 93, 97, 114};
	knownDFill[93] = {
		93, 2, 92, 9, 87, 12, 84, 31, 76, 35, 64, 39, 53, 94, 98, 115};
	knownDFill[94] = {94, 2, 93, 9, 88, 13, 84, 31, 77, 35, 95, 99, 116};
	knownDFill[95] = {95, 2, 94, 9, 89, 14, 84, 31, 35, 117};
	knownDFill[96] = {
		96, 2, 9, 91, 16, 84, 31, 78, 35, 65, 39, 54, 45, 46, 97, 100, 118};
	knownDFill[97] = {
		97, 2, 96, 9, 92, 17, 84, 31, 79, 35, 66, 39, 55, 98, 101, 119};
	knownDFill[98] = {98, 2, 97, 9, 93, 18, 84, 31, 80, 35, 99, 102, 120};
	knownDFill[99] = {99, 2, 98, 9, 94, 19, 84, 31, 35, 121};
	knownDFill[100] = {
		100, 2, 9, 96, 21, 84, 31, 81, 35, 67, 39, 56, 101, 103, 122};
	knownDFill[101] = {101, 2, 100, 9, 97, 22, 84, 31, 82, 35, 102, 104, 123};
	knownDFill[102] = {102, 2, 101, 9, 98, 23, 84, 31, 35, 124};
	knownDFill[103] = {103, 2, 9, 100, 25, 84, 31, 83, 35, 104, 105, 125};
	knownDFill[104] = {104, 2, 103, 9, 101, 26, 84, 31, 35, 126};
	knownDFill[105] = {105, 9, 103, 28, 84, 31, 35, 127};
	knownDFill[106] = {106, 0, 107, 1, 31, 84, 35, 68, 39, 57, 46};
	knownDFill[107] = {107, 0, 1, 106, 2, 9, 31, 85, 32, 84, 35, 69, 36, 68, 39,
		58, 40, 57, 46, 47, 108, 113, 128};
	knownDFill[108] = {108, 2, 107, 9, 31, 86, 33, 84, 35, 70, 37, 68, 39, 59,
		41, 57, 46, 48, 109, 114, 129};
	knownDFill[109] = {
		109, 2, 108, 9, 31, 87, 35, 71, 39, 60, 42, 57, 46, 49, 110, 115, 130};
	knownDFill[110] = {
		110, 2, 109, 9, 31, 88, 35, 72, 39, 61, 46, 50, 111, 116, 131};
	knownDFill[111] = {111, 2, 110, 9, 31, 89, 35, 73, 112, 117, 132};
	knownDFill[112] = {112, 2, 111, 31, 90, 133};
	knownDFill[113] = {113, 2, 9, 107, 31, 91, 34, 84, 35, 74, 38, 68, 39, 62,
		43, 57, 46, 51, 114, 118, 135};
	knownDFill[114] = {114, 2, 113, 9, 108, 31, 92, 35, 75, 39, 63, 44, 57, 46,
		52, 115, 119, 136};
	knownDFill[115] = {
		115, 2, 114, 9, 109, 31, 93, 35, 76, 39, 64, 46, 53, 116, 120, 137};
	knownDFill[116] = {116, 2, 115, 9, 110, 31, 94, 35, 77, 117, 121, 138};
	knownDFill[117] = {117, 2, 116, 9, 111, 31, 95, 139};
	knownDFill[118] = {
		118, 2, 9, 113, 31, 96, 35, 78, 39, 65, 45, 57, 46, 54, 119, 122, 141};
	knownDFill[119] = {
		119, 2, 118, 9, 114, 31, 97, 35, 79, 39, 66, 46, 55, 120, 123, 142};
	knownDFill[120] = {120, 2, 119, 9, 115, 31, 98, 35, 80, 121, 124, 143};
	knownDFill[121] = {121, 2, 120, 9, 116, 31, 99, 144};
	knownDFill[122] = {
		122, 2, 9, 118, 31, 100, 35, 81, 39, 67, 46, 56, 123, 125, 146};
	knownDFill[123] = {123, 2, 122, 9, 119, 31, 101, 35, 82, 124, 126, 147};
	knownDFill[124] = {124, 2, 123, 9, 120, 31, 102, 148};
	knownDFill[125] = {125, 2, 9, 122, 31, 103, 35, 83, 126, 127, 150};
	knownDFill[126] = {126, 2, 125, 9, 123, 31, 104, 151};
	knownDFill[127] = {127, 9, 125, 31, 105, 153};
	knownDFill[128] = {128, 2, 9, 31, 106, 107, 35, 84, 85, 36, 39, 68, 69, 40,
		46, 57, 58, 47, 129, 135};
	knownDFill[129] = {129, 2, 128, 9, 31, 108, 35, 86, 37, 84, 39, 70, 41, 68,
		46, 59, 48, 57, 130, 136};
	knownDFill[130] = {130, 2, 129, 9, 31, 109, 35, 87, 39, 71, 42, 68, 46, 60,
		49, 57, 131, 137};
	knownDFill[131] = {
		131, 2, 130, 9, 31, 110, 35, 88, 39, 72, 46, 61, 50, 57, 132, 138};
	knownDFill[132] = {132, 2, 131, 9, 31, 111, 35, 89, 39, 73, 133, 139};
	knownDFill[133] = {133, 2, 132, 9, 31, 112, 35, 90, 134, 140};
	knownDFill[134] = {134, 2, 133};
	knownDFill[135] = {135, 2, 9, 128, 31, 113, 35, 91, 38, 84, 39, 74, 43, 68,
		46, 62, 51, 57, 136, 141};
	knownDFill[136] = {136, 2, 135, 9, 129, 31, 114, 35, 92, 39, 75, 44, 68, 46,
		63, 52, 57, 137, 142};
	knownDFill[137] = {
		137, 2, 136, 9, 130, 31, 115, 35, 93, 39, 76, 46, 64, 53, 57, 138, 143};
	knownDFill[138] = {138, 2, 137, 9, 131, 31, 116, 35, 94, 39, 77, 139, 144};
	knownDFill[139] = {139, 2, 138, 9, 132, 31, 117, 35, 95, 140, 145};
	knownDFill[140] = {140, 2, 139, 9, 133};
	knownDFill[141] = {141, 2, 9, 135, 31, 118, 35, 96, 39, 78, 45, 68, 46, 65,
		54, 57, 142, 146};
	knownDFill[142] = {
		142, 2, 141, 9, 136, 31, 119, 35, 97, 39, 79, 46, 66, 55, 57, 143, 147};
	knownDFill[143] = {143, 2, 142, 9, 137, 31, 120, 35, 98, 39, 80, 144, 148};
	knownDFill[144] = {144, 2, 143, 9, 138, 31, 121, 35, 99, 145, 149};
	knownDFill[145] = {145, 2, 144, 9, 139};
	knownDFill[146] = {
		146, 2, 9, 141, 31, 122, 35, 100, 39, 81, 46, 67, 56, 57, 147, 150};
	knownDFill[147] = {147, 2, 146, 9, 142, 31, 123, 35, 101, 39, 82, 148, 151};
	knownDFill[148] = {148, 2, 147, 9, 143, 31, 124, 35, 102, 149, 152};
	knownDFill[149] = {149, 2, 148, 9, 144};
	knownDFill[150] = {150, 2, 9, 146, 31, 125, 35, 103, 39, 83, 151, 153};
	knownDFill[151] = {151, 2, 150, 9, 147, 31, 126, 35, 104, 152, 154};
	knownDFill[152] = {152, 2, 151, 9, 148};
	knownDFill[153] = {153, 2, 9, 150, 31, 127, 35, 105, 154, 155};
	knownDFill[154] = {154, 2, 153, 9, 151};
	knownDFill[155] = {155, 9, 153};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 2547);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 290600199485, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::He, 2), 121.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::D, 1), 229.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::I, 1), 0.15785, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::T, 1), 228.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 3196583454390, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(
			fluxes[i], xolotl::test::fullyRefined::knownFluxes[i], 0.01);
	}

	// Check the partials computation
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	int startingIdx = 0;
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
				auto index = std::distance(row.begin(), iter);
				BOOST_REQUIRE_CLOSE(hPartials[startingIdx + index],
					xolotl::test::fullyRefined::knownPartials[startingIdx + j],
					0.01);
			}
			startingIdx += row.size();
		}
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::V] = 1;
	auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 1);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.137265, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 3.6, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 505312.69, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 1800000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.3, 0.01);
	const auto& reg0 = cluster.getRegion();
	Composition lo = reg0.getOrigin();
	Composition hi = reg0.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 4);

	comp[Spec::V] = 0;
	comp[Spec::I] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.15785, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 10.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 78358278338, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 88000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.01, 0.01);
	const auto& reg1 = cluster.getRegion();
	lo = reg1.getOrigin();
	hi = reg1.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 4);

	comp[Spec::I] = 0;
	comp[Spec::He] = 5;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 57);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3648, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 26.1, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 1242214406, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 5000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.12, 0.01);
	const auto& reg2 = cluster.getRegion();
	lo = reg2.getOrigin();
	hi = reg2.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 5);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 6);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 4);

	comp[Spec::He] = 0;
	comp[Spec::D] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 9);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.075, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getFormationEnergy(), std::numeric_limits<double>::infinity());
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 3440887974, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 283000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.38, 0.01);
	const auto& reg3 = cluster.getRegion();
	lo = reg3.getOrigin();
	hi = reg3.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 4);

	comp[Spec::D] = 0;
	comp[Spec::T] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 2);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.075, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getFormationEnergy(), std::numeric_limits<double>::infinity());
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 2808640007, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 231000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.38, 0.01);
	const auto& reg4 = cluster.getRegion();
	lo = reg4.getOrigin();
	hi = reg4.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 2);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 4);

	comp[Spec::V] = 1;
	comp[Spec::He] = 8;
	comp[Spec::D] = 3;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 123);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.137265, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 30.1049, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg5 = cluster.getRegion();
	lo = reg5.getOrigin();
	hi = reg5.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 8);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 9);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 3);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 4);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 2);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 4);
}

BOOST_AUTO_TEST_CASE(reducedMatrixMethod)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 1 1 1 1" << std::endl
			  << "process=reaction" << std::endl
			  << "tempParam=1000" << std::endl
			  << "petscArgs=-snes_mf_operator" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIFullSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the temperature and lattice constant from the options
	double latticeConst = opts.getLatticeParameter() <= 0.0 ?
		xolotl::core::tungstenLatticeConstant :
		opts.getLatticeParameter();
	double temperature = opts.getTempParam();

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe =
		util::getMaxHePerVLoop(maxV, latticeConst, temperature);
	NetworkType::AmountType maxD = 2.0 / 3.0 * (double)maxHe;
	NetworkType::AmountType maxT = 2.0 / 3.0 * (double)maxHe;
	NetworkType network({maxHe, maxD, maxT, maxV, maxI}, 1, opts);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0};
	knownDFill[1] = {1};
	knownDFill[2] = {2};
	knownDFill[3] = {3};
	knownDFill[4] = {4};
	knownDFill[5] = {5};
	knownDFill[6] = {6};
	knownDFill[7] = {7};
	knownDFill[8] = {8};
	knownDFill[9] = {9};
	knownDFill[10] = {10};
	knownDFill[11] = {11};
	knownDFill[12] = {12};
	knownDFill[13] = {13};
	knownDFill[14] = {14};
	knownDFill[15] = {15};
	knownDFill[16] = {16};
	knownDFill[17] = {17};
	knownDFill[18] = {18};
	knownDFill[19] = {19};
	knownDFill[20] = {20};
	knownDFill[21] = {21};
	knownDFill[22] = {22};
	knownDFill[23] = {23};
	knownDFill[24] = {24};
	knownDFill[25] = {25};
	knownDFill[26] = {26};
	knownDFill[27] = {27};
	knownDFill[28] = {28};
	knownDFill[29] = {29};
	knownDFill[30] = {30};
	knownDFill[31] = {31};
	knownDFill[32] = {32};
	knownDFill[33] = {33};
	knownDFill[34] = {34};
	knownDFill[35] = {35};
	knownDFill[36] = {36};
	knownDFill[37] = {37};
	knownDFill[38] = {38};
	knownDFill[39] = {39};
	knownDFill[40] = {40};
	knownDFill[41] = {41};
	knownDFill[42] = {42};
	knownDFill[43] = {43};
	knownDFill[44] = {44};
	knownDFill[45] = {45};
	knownDFill[46] = {46};
	knownDFill[47] = {47};
	knownDFill[48] = {48};
	knownDFill[49] = {49};
	knownDFill[50] = {50};
	knownDFill[51] = {51};
	knownDFill[52] = {52};
	knownDFill[53] = {53};
	knownDFill[54] = {54};
	knownDFill[55] = {55};
	knownDFill[56] = {56};
	knownDFill[57] = {57};
	knownDFill[58] = {58};
	knownDFill[59] = {59};
	knownDFill[60] = {60};
	knownDFill[61] = {61};
	knownDFill[62] = {62};
	knownDFill[63] = {63};
	knownDFill[64] = {64};
	knownDFill[65] = {65};
	knownDFill[66] = {66};
	knownDFill[67] = {67};
	knownDFill[68] = {68};
	knownDFill[69] = {69};
	knownDFill[70] = {70};
	knownDFill[71] = {71};
	knownDFill[72] = {72};
	knownDFill[73] = {73};
	knownDFill[74] = {74};
	knownDFill[75] = {75};
	knownDFill[76] = {76};
	knownDFill[77] = {77};
	knownDFill[78] = {78};
	knownDFill[79] = {79};
	knownDFill[80] = {80};
	knownDFill[81] = {81};
	knownDFill[82] = {82};
	knownDFill[83] = {83};
	knownDFill[84] = {84};
	knownDFill[85] = {85};
	knownDFill[86] = {86};
	knownDFill[87] = {87};
	knownDFill[88] = {88};
	knownDFill[89] = {89};
	knownDFill[90] = {90};
	knownDFill[91] = {91};
	knownDFill[92] = {92};
	knownDFill[93] = {93};
	knownDFill[94] = {94};
	knownDFill[95] = {95};
	knownDFill[96] = {96};
	knownDFill[97] = {97};
	knownDFill[98] = {98};
	knownDFill[99] = {99};
	knownDFill[100] = {100};
	knownDFill[101] = {101};
	knownDFill[102] = {102};
	knownDFill[103] = {103};
	knownDFill[104] = {104};
	knownDFill[105] = {105};
	knownDFill[106] = {106};
	knownDFill[107] = {107};
	knownDFill[108] = {108};
	knownDFill[109] = {109};
	knownDFill[110] = {110};
	knownDFill[111] = {111};
	knownDFill[112] = {112};
	knownDFill[113] = {113};
	knownDFill[114] = {114};
	knownDFill[115] = {115};
	knownDFill[116] = {116};
	knownDFill[117] = {117};
	knownDFill[118] = {118};
	knownDFill[119] = {119};
	knownDFill[120] = {120};
	knownDFill[121] = {121};
	knownDFill[122] = {122};
	knownDFill[123] = {123};
	knownDFill[124] = {124};
	knownDFill[125] = {125};
	knownDFill[126] = {126};
	knownDFill[127] = {127};
	knownDFill[128] = {128};
	knownDFill[129] = {129};
	knownDFill[130] = {130};
	knownDFill[131] = {131};
	knownDFill[132] = {132};
	knownDFill[133] = {133};
	knownDFill[134] = {134};
	knownDFill[135] = {135};
	knownDFill[136] = {136};
	knownDFill[137] = {137};
	knownDFill[138] = {138};
	knownDFill[139] = {139};
	knownDFill[140] = {140};
	knownDFill[141] = {141};
	knownDFill[142] = {142};
	knownDFill[143] = {143};
	knownDFill[144] = {144};
	knownDFill[145] = {145};
	knownDFill[146] = {146};
	knownDFill[147] = {147};
	knownDFill[148] = {148};
	knownDFill[149] = {149};
	knownDFill[150] = {150};
	knownDFill[151] = {151};
	knownDFill[152] = {152};
	knownDFill[153] = {153};
	knownDFill[154] = {154};
	knownDFill[155] = {155};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 156);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 290600199485, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::He, 2), 121.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::D, 1), 229.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalRadiusConcentration(dConcs, Spec::I, 1), 0.15785, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::T, 1), 228.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 3196583454390, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(
			fluxes[i], xolotl::test::reducedMatrixMethod::knownFluxes[i], 0.01);
	}

	// Check the partials computation
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	int startingIdx = 0;
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
				auto index = std::distance(row.begin(), iter);
				BOOST_REQUIRE_CLOSE(hPartials[startingIdx + index],
					xolotl::test::reducedMatrixMethod::knownPartials
						[startingIdx + j],
					0.01);
			}
			startingIdx += row.size();
		}
	}
}

BOOST_AUTO_TEST_CASE(HeliumSpeciesList)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 0 2 2" << std::endl
			  << "process=reaction" << std::endl
			  << "tempParam=1000" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIHeliumSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the temperature and lattice constant from the options
	double latticeConst = opts.getLatticeParameter() <= 0.0 ?
		xolotl::core::tungstenLatticeConstant :
		opts.getLatticeParameter();
	double temperature = opts.getTempParam();

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe =
		util::getMaxHePerVLoop(maxV, latticeConst, temperature);
	NetworkType network({maxHe, maxV, maxI}, 1, opts);

	BOOST_REQUIRE(!network.hasDeuterium());
	BOOST_REQUIRE(!network.hasTritium());

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 35);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 35);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.317, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.0159275, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 35);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 3);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 3);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 2);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 2, 3, 5, 6, 8, 9, 11, 12, 14, 15, 17, 18, 20, 21, 23,
		24, 26, 27, 29, 1, 4, 25, 28, 7, 22, 10, 19, 13, 16};
	knownDFill[1] = {
		1, 0, 2, 3, 6, 9, 12, 15, 18, 21, 24, 27, 7, 25, 10, 22, 13, 19, 16};
	knownDFill[2] = {2, 0, 3, 1, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22,
		23, 25, 26, 28, 6, 9, 12, 15, 18, 21, 24, 27, 29};
	knownDFill[3] = {3, 0, 1, 2, 4, 7, 10, 13, 16, 19, 22, 6};
	knownDFill[4] = {4, 0, 5, 1, 6, 2, 3, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
		17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
	knownDFill[5] = {5, 0, 6, 2, 4, 7, 10, 13, 16, 19, 22, 8};
	knownDFill[6] = {6, 0, 1, 2, 5, 3, 4, 7, 10, 13, 16, 19, 22, 9};
	knownDFill[7] = {7, 0, 8, 1, 9, 2, 3, 4, 5, 6, 10, 11, 12, 13, 14, 15, 16,
		17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32};
	knownDFill[8] = {8, 0, 9, 2, 7, 4, 5, 10, 13, 16, 19, 22, 11};
	knownDFill[9] = {9, 0, 1, 2, 8, 3, 7, 4, 6, 10, 13, 16, 19, 22, 12};
	knownDFill[10] = {10, 0, 11, 1, 12, 2, 3, 4, 7, 5, 6, 8, 9, 13, 14, 15, 16,
		17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31};
	knownDFill[11] = {11, 0, 12, 2, 10, 4, 8, 5, 7, 13, 16, 19, 22, 14};
	knownDFill[12] = {12, 0, 1, 2, 11, 3, 10, 4, 9, 6, 7, 13, 16, 19, 22, 15};
	knownDFill[13] = {13, 0, 14, 1, 15, 2, 3, 4, 10, 5, 6, 7, 8, 9, 11, 12, 16,
		17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
	knownDFill[14] = {14, 0, 15, 2, 13, 4, 11, 5, 10, 7, 8, 16, 19, 22, 17};
	knownDFill[15] = {
		15, 0, 1, 2, 14, 3, 13, 4, 12, 6, 10, 7, 9, 16, 19, 22, 18};
	knownDFill[16] = {16, 0, 17, 1, 18, 2, 3, 4, 13, 5, 6, 7, 10, 8, 9, 11, 12,
		14, 15, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29};
	knownDFill[17] = {17, 0, 18, 2, 16, 4, 14, 5, 13, 7, 11, 8, 10, 19, 22, 20};
	knownDFill[18] = {
		18, 0, 1, 2, 17, 3, 16, 4, 15, 6, 13, 7, 12, 9, 10, 19, 22, 21};
	knownDFill[19] = {19, 0, 20, 1, 21, 2, 3, 4, 16, 5, 6, 7, 13, 8, 9, 10, 11,
		12, 14, 15, 17, 18, 22, 23, 24, 25, 26, 27};
	knownDFill[20] = {
		20, 0, 21, 2, 19, 4, 17, 5, 16, 7, 14, 8, 13, 10, 11, 22, 23};
	knownDFill[21] = {
		21, 0, 1, 2, 20, 3, 19, 4, 18, 6, 16, 7, 15, 9, 13, 10, 12, 22, 24};
	knownDFill[22] = {22, 0, 23, 1, 24, 2, 3, 4, 19, 5, 6, 7, 16, 8, 9, 10, 13,
		11, 12, 14, 15, 17, 18, 20, 21, 25};
	knownDFill[23] = {
		23, 0, 24, 2, 22, 4, 20, 5, 19, 7, 17, 8, 16, 10, 14, 11, 13, 26};
	knownDFill[24] = {
		24, 0, 1, 2, 23, 3, 22, 4, 21, 6, 19, 7, 18, 9, 16, 10, 15, 12, 13, 27};
	knownDFill[25] = {25, 0, 26, 1, 27, 2, 4, 22, 7, 19, 10, 16, 13};
	knownDFill[26] = {26, 0, 27, 2, 25, 4, 23, 5, 22, 7, 20, 8, 19, 10, 17, 11,
		16, 13, 14, 28};
	knownDFill[27] = {27, 0, 1, 2, 26, 4, 24, 6, 22, 7, 21, 9, 19, 10, 18, 12,
		16, 13, 15, 29};
	knownDFill[28] = {
		28, 0, 29, 2, 4, 25, 26, 7, 22, 23, 8, 10, 19, 20, 11, 13, 16, 17, 14};
	knownDFill[29] = {
		29, 0, 2, 28, 4, 27, 7, 24, 9, 22, 10, 21, 12, 19, 13, 18, 15, 16, 30};
	knownDFill[30] = {30, 4, 28, 29, 7, 25, 26, 27, 10, 22, 23, 24, 11, 12, 13,
		19, 20, 21, 14, 15, 16, 17, 18, 31};
	knownDFill[31] = {31, 4, 30, 7, 28, 29, 10, 25, 26, 27, 13, 22, 23, 24, 14,
		15, 16, 19, 20, 21, 17, 18, 32};
	knownDFill[32] = {32, 4, 31, 7, 30, 10, 28, 29, 13, 25, 26, 27, 16, 22, 23,
		24, 17, 18, 19, 20, 21, 33};
	knownDFill[33] = {33, 4, 32, 7, 31, 10, 30, 13, 28, 29, 16, 25, 26, 27, 19,
		22, 23, 24, 20, 21, 34};
	knownDFill[34] = {34, 4, 33, 7, 32, 10, 31, 13, 30, 16, 28, 29, 19, 25, 26,
		27, 22, 23, 24};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 730);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 621745254537, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::He, 2), 28.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::V, 1), 40.0, 0.01);
	BOOST_REQUIRE_CLOSE(network.getTotalRadiusConcentration(dConcs, Spec::I, 1),
		0.356278, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::He, 1), 150.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 6542589997348, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(
			fluxes[i], xolotl::test::HeliumSpeciesList::knownFluxes[i], 0.01);
	}

	// Check the partials computation
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	int startingIdx = 0;
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
				auto index = std::distance(row.begin(), iter);
				BOOST_REQUIRE_CLOSE(hPartials[startingIdx + index],
					xolotl::test::HeliumSpeciesList::knownPartials[startingIdx +
						j],
					0.01);
			}
			startingIdx += row.size();
		}
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::V] = 1;
	auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 2);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.137265, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 3.6, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 505312.69, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 1800000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.3, 0.01);
	const auto& reg0 = cluster.getRegion();
	Composition lo = reg0.getOrigin();
	Composition hi = reg0.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::V] = 0;
	comp[Spec::I] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.15785, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 10.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 78358278338, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 88000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.01, 0.01);
	const auto& reg1 = cluster.getRegion();
	lo = reg1.getOrigin();
	hi = reg1.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::I] = 0;
	comp[Spec::He] = 5;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 16);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3648, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 26.1, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 1242214406, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 5000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.12, 0.01);
	const auto& reg2 = cluster.getRegion();
	lo = reg2.getOrigin();
	hi = reg2.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 5);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 6);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);

	comp[Spec::He] = 4;
	comp[Spec::V] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 14);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.137265, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 14.8829, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg3 = cluster.getRegion();
	lo = reg3.getOrigin();
	hi = reg3.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 4);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 5);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 2);
}

BOOST_AUTO_TEST_CASE(DeuteriumSpeciesList)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 1 0 1 1" << std::endl
			  << "process=reaction" << std::endl
			  << "tempParam=1000" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIDeuteriumSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the temperature and lattice constant from the options
	double latticeConst = opts.getLatticeParameter() <= 0.0 ?
		xolotl::core::tungstenLatticeConstant :
		opts.getLatticeParameter();
	double temperature = opts.getTempParam();

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe =
		util::getMaxHePerVLoop(maxV, latticeConst, temperature);
	NetworkType::AmountType maxD = 2.0 / 3.0 * (double)maxHe;
	NetworkType network({maxHe, maxD, maxV, maxI}, 1, opts);

	BOOST_REQUIRE(network.hasDeuterium());
	BOOST_REQUIRE(!network.hasTritium());

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 56);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 56);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.317, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.0159275, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 56);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 4);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 4);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 3);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {
		0, 1, 3, 10, 13, 16, 20, 25, 30, 36, 43, 9, 42, 12, 35, 15, 29, 19, 24};
	knownDFill[1] = {1, 0, 2, 9, 12, 15, 19, 24, 29, 35, 42, 3, 10, 13, 16, 20,
		25, 30, 36, 43};
	knownDFill[2] = {2, 0, 3, 1, 4, 5, 6, 7, 10, 13, 16, 17, 20, 21, 22, 25, 26,
		27, 30, 31, 32, 33, 36, 37, 38, 39, 40, 43, 44, 45, 46, 47, 49, 50, 51,
		52, 53, 54, 8, 11, 14, 18, 23, 28, 34, 41, 48, 55};
	knownDFill[3] = {3, 0, 1, 2, 9, 12, 15, 19, 24, 29, 35, 4, 11};
	knownDFill[4] = {4, 2, 3, 15, 19, 24, 29, 35, 5};
	knownDFill[5] = {5, 2, 4, 19, 24, 29, 35, 6};
	knownDFill[6] = {6, 2, 5, 29, 35, 7};
	knownDFill[7] = {7, 2, 6, 35, 8};
	knownDFill[8] = {8, 2, 7};
	knownDFill[9] = {9, 0, 10, 1, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
		22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
		40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54};
	knownDFill[10] = {10, 0, 1, 9, 2, 12, 15, 19, 24, 29, 35, 11, 13};
	knownDFill[11] = {11, 2, 10, 3, 9, 12, 15, 19, 24, 29, 35, 14};
	knownDFill[12] = {12, 0, 13, 1, 3, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20,
		21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
		39, 40, 41};
	knownDFill[13] = {13, 0, 1, 12, 2, 9, 10, 15, 19, 24, 29, 35, 14, 16};
	knownDFill[14] = {14, 2, 13, 3, 12, 9, 11, 15, 19, 24, 29, 35, 17};
	knownDFill[15] = {15, 0, 16, 1, 3, 4, 9, 12, 10, 11, 13, 14, 17, 18, 19, 20,
		21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
	knownDFill[16] = {16, 0, 1, 15, 2, 9, 13, 10, 12, 19, 24, 29, 17, 20};
	knownDFill[17] = {17, 2, 16, 3, 15, 9, 14, 11, 12, 19, 24, 29, 18, 21};
	knownDFill[18] = {18, 2, 17, 4, 15, 9, 12, 19, 24, 29, 22};
	knownDFill[19] = {19, 0, 20, 1, 3, 4, 5, 9, 15, 10, 11, 12, 13, 14, 16, 17,
		18, 21, 22, 23, 24, 25, 26, 27, 28};
	knownDFill[20] = {20, 0, 1, 19, 2, 9, 16, 10, 15, 12, 13, 24, 21, 25};
	knownDFill[21] = {21, 2, 20, 3, 19, 9, 17, 11, 15, 12, 14, 24, 22, 26};
	knownDFill[22] = {22, 2, 21, 4, 19, 9, 18, 12, 15, 24, 23, 27};
	knownDFill[23] = {23, 2, 22, 5, 19, 9, 12, 15, 24, 28};
	knownDFill[24] = {24, 0, 25, 1, 3, 4, 5, 9, 19, 10, 11, 12, 15, 13, 14, 16,
		17, 18, 20, 21, 22, 23, 29};
	knownDFill[25] = {25, 0, 1, 24, 2, 9, 20, 10, 19, 12, 16, 13, 15, 26, 30};
	knownDFill[26] = {26, 2, 25, 3, 24, 9, 21, 11, 19, 12, 17, 14, 15, 27, 31};
	knownDFill[27] = {27, 2, 26, 4, 24, 9, 22, 12, 18, 15, 19, 28, 32};
	knownDFill[28] = {28, 2, 27, 5, 24, 9, 23, 12, 15, 19, 33};
	knownDFill[29] = {29, 0, 30, 1, 3, 4, 5, 6, 9, 24, 10, 11, 12, 19, 13, 14,
		15, 16, 17, 18, 35};
	knownDFill[30] = {
		30, 0, 1, 29, 2, 9, 25, 10, 24, 12, 20, 13, 19, 15, 16, 31, 36};
	knownDFill[31] = {
		31, 2, 30, 3, 29, 9, 26, 11, 24, 12, 21, 14, 19, 15, 17, 32, 37};
	knownDFill[32] = {32, 2, 31, 4, 29, 9, 27, 12, 22, 15, 18, 33, 38};
	knownDFill[33] = {33, 2, 32, 5, 29, 9, 28, 12, 23, 15, 34, 39};
	knownDFill[34] = {34, 2, 33, 6, 29, 9, 12, 15, 40};
	knownDFill[35] = {
		35, 0, 36, 1, 3, 4, 5, 6, 7, 9, 29, 10, 11, 12, 24, 13, 14, 15, 19, 42};
	knownDFill[36] = {
		36, 0, 1, 35, 2, 9, 30, 10, 29, 12, 25, 13, 24, 15, 20, 16, 19, 37, 43};
	knownDFill[37] = {37, 2, 36, 3, 35, 9, 31, 11, 29, 12, 26, 14, 24, 15, 21,
		17, 19, 38, 44};
	knownDFill[38] = {38, 2, 37, 4, 35, 9, 32, 12, 27, 15, 22, 18, 19, 39, 45};
	knownDFill[39] = {39, 2, 38, 5, 35, 9, 33, 12, 28, 15, 23, 40, 46};
	knownDFill[40] = {40, 2, 39, 6, 35, 9, 34, 12, 41, 47};
	knownDFill[41] = {41, 2, 40, 7, 35, 9, 12, 48};
	knownDFill[42] = {42, 0, 43, 1, 9, 35, 12, 29, 15, 24, 19};
	knownDFill[43] = {43, 0, 1, 42, 2, 9, 36, 10, 35, 12, 30, 13, 29, 15, 25,
		16, 24, 19, 20, 44, 49};
	knownDFill[44] = {44, 2, 43, 9, 37, 11, 35, 12, 31, 14, 29, 15, 26, 17, 24,
		19, 21, 45, 50};
	knownDFill[45] = {45, 2, 44, 9, 38, 12, 32, 15, 27, 18, 24, 19, 22, 46, 51};
	knownDFill[46] = {46, 2, 45, 9, 39, 12, 33, 15, 28, 19, 23, 47, 52};
	knownDFill[47] = {47, 2, 46, 9, 40, 12, 34, 48, 53};
	knownDFill[48] = {48, 2, 47, 9, 41, 54};
	knownDFill[49] = {
		49, 2, 9, 42, 43, 12, 35, 36, 13, 15, 29, 30, 16, 19, 24, 25, 20, 50};
	knownDFill[50] = {
		50, 2, 49, 9, 44, 12, 37, 14, 35, 15, 31, 17, 29, 19, 26, 21, 24, 51};
	knownDFill[51] = {
		51, 2, 50, 9, 45, 12, 38, 15, 32, 18, 29, 19, 27, 22, 24, 52};
	knownDFill[52] = {52, 2, 51, 9, 46, 12, 39, 15, 33, 19, 28, 23, 24, 53};
	knownDFill[53] = {53, 2, 52, 9, 47, 12, 40, 15, 34, 54};
	knownDFill[54] = {54, 2, 53, 9, 48, 12, 41, 55};
	knownDFill[55] = {55, 2, 54};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 863);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 290600199485, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::He, 2), 44.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::V, 1), 46.0, 0.01);
	BOOST_REQUIRE_CLOSE(network.getTotalRadiusConcentration(dConcs, Spec::I, 1),
		0.15785478, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::D, 1), 99.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 2905985128900, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(fluxes[i],
			xolotl::test::DeuteriumSpeciesList::knownFluxes[i], 0.01);
	}

	// Check the partials computation
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	int startingIdx = 0;
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
				auto index = std::distance(row.begin(), iter);
				BOOST_REQUIRE_CLOSE(hPartials[startingIdx + index],
					xolotl::test::DeuteriumSpeciesList::knownPartials
						[startingIdx + j],
					0.01);
			}
			startingIdx += row.size();
		}
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::V] = 1;
	auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 1);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.137265, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 3.6, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 505312.69, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 1800000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.3, 0.01);
	const auto& reg0 = cluster.getRegion();
	Composition lo = reg0.getOrigin();
	Composition hi = reg0.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);

	comp[Spec::V] = 0;
	comp[Spec::I] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.15785, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 10.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 78358278338, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 88000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.01, 0.01);
	const auto& reg1 = cluster.getRegion();
	lo = reg1.getOrigin();
	hi = reg1.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);

	comp[Spec::I] = 0;
	comp[Spec::He] = 5;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 24);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3648, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 26.1, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 1242214406, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 5000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.12, 0.01);
	const auto& reg2 = cluster.getRegion();
	lo = reg2.getOrigin();
	hi = reg2.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 5);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 6);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);

	comp[Spec::He] = 0;
	comp[Spec::D] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 2);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.075, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getFormationEnergy(), std::numeric_limits<double>::infinity());
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 3440887974, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 283000000000, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.38, 0.01);
	const auto& reg3 = cluster.getRegion();
	lo = reg3.getOrigin();
	hi = reg3.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 2);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);

	comp[Spec::He] = 3;
	comp[Spec::V] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 17);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.137265, 0.01);
	BOOST_REQUIRE_EQUAL(cluster.getFormationEnergy(), 11.5304);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg4 = cluster.getRegion();
	lo = reg4.getOrigin();
	hi = reg4.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 3);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 4);
	BOOST_REQUIRE_EQUAL(lo[Spec::D], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::D], 2);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);
}

BOOST_AUTO_TEST_CASE(TritiumSpeciesList)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 1 1 1" << std::endl
			  << "process=reaction" << std::endl
			  << "tempParam=1000" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSITritiumSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the temperature and lattice constant from the options
	double latticeConst = opts.getLatticeParameter() <= 0.0 ?
		xolotl::core::tungstenLatticeConstant :
		opts.getLatticeParameter();
	double temperature = opts.getTempParam();

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe =
		util::getMaxHePerVLoop(maxV, latticeConst, temperature);
	NetworkType::AmountType maxT = 2.0 / 3.0 * (double)maxHe;
	NetworkType network({maxHe, maxT, maxV, maxI}, 1, opts);

	BOOST_REQUIRE(!network.hasDeuterium());
	BOOST_REQUIRE(network.hasTritium());

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 56);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 56);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.317, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.0159275, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableReSolution() == false);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 56);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 4);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 4);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 3);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {
		0, 1, 3, 10, 13, 16, 20, 25, 30, 36, 43, 9, 42, 12, 35, 15, 29, 19, 24};
	knownDFill[1] = {1, 0, 2, 9, 12, 15, 19, 24, 29, 35, 42, 3, 10, 13, 16, 20,
		25, 30, 36, 43};
	knownDFill[2] = {2, 0, 3, 1, 4, 5, 6, 7, 10, 13, 16, 17, 20, 21, 22, 25, 26,
		27, 30, 31, 32, 33, 36, 37, 38, 39, 40, 43, 44, 45, 46, 47, 49, 50, 51,
		52, 53, 54, 8, 11, 14, 18, 23, 28, 34, 41, 48, 55};
	knownDFill[3] = {3, 0, 1, 2, 9, 12, 15, 19, 24, 29, 35, 4, 11};
	knownDFill[4] = {4, 2, 3, 15, 19, 24, 29, 35, 5};
	knownDFill[5] = {5, 2, 4, 19, 24, 29, 35, 6};
	knownDFill[6] = {6, 2, 5, 29, 35, 7};
	knownDFill[7] = {7, 2, 6, 35, 8};
	knownDFill[8] = {8, 2, 7};
	knownDFill[9] = {9, 0, 10, 1, 3, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
		22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,
		40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54};
	knownDFill[10] = {10, 0, 1, 9, 2, 12, 15, 19, 24, 29, 35, 11, 13};
	knownDFill[11] = {11, 2, 10, 3, 9, 12, 15, 19, 24, 29, 35, 14};
	knownDFill[12] = {12, 0, 13, 1, 3, 9, 10, 11, 14, 15, 16, 17, 18, 19, 20,
		21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
		39, 40, 41};
	knownDFill[13] = {13, 0, 1, 12, 2, 9, 10, 15, 19, 24, 29, 35, 14, 16};
	knownDFill[14] = {14, 2, 13, 3, 12, 9, 11, 15, 19, 24, 29, 35, 17};
	knownDFill[15] = {15, 0, 16, 1, 3, 4, 9, 12, 10, 11, 13, 14, 17, 18, 19, 20,
		21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
	knownDFill[16] = {16, 0, 1, 15, 2, 9, 13, 10, 12, 19, 24, 29, 17, 20};
	knownDFill[17] = {17, 2, 16, 3, 15, 9, 14, 11, 12, 19, 24, 29, 18, 21};
	knownDFill[18] = {18, 2, 17, 4, 15, 9, 12, 19, 24, 29, 22};
	knownDFill[19] = {19, 0, 20, 1, 3, 4, 5, 9, 15, 10, 11, 12, 13, 14, 16, 17,
		18, 21, 22, 23, 24, 25, 26, 27, 28};
	knownDFill[20] = {20, 0, 1, 19, 2, 9, 16, 10, 15, 12, 13, 24, 21, 25};
	knownDFill[21] = {21, 2, 20, 3, 19, 9, 17, 11, 15, 12, 14, 24, 22, 26};
	knownDFill[22] = {22, 2, 21, 4, 19, 9, 18, 12, 15, 24, 23, 27};
	knownDFill[23] = {23, 2, 22, 5, 19, 9, 12, 15, 24, 28};
	knownDFill[24] = {24, 0, 25, 1, 3, 4, 5, 9, 19, 10, 11, 12, 15, 13, 14, 16,
		17, 18, 20, 21, 22, 23, 29};
	knownDFill[25] = {25, 0, 1, 24, 2, 9, 20, 10, 19, 12, 16, 13, 15, 26, 30};
	knownDFill[26] = {26, 2, 25, 3, 24, 9, 21, 11, 19, 12, 17, 14, 15, 27, 31};
	knownDFill[27] = {27, 2, 26, 4, 24, 9, 22, 12, 18, 15, 19, 28, 32};
	knownDFill[28] = {28, 2, 27, 5, 24, 9, 23, 12, 15, 19, 33};
	knownDFill[29] = {29, 0, 30, 1, 3, 4, 5, 6, 9, 24, 10, 11, 12, 19, 13, 14,
		15, 16, 17, 18, 35};
	knownDFill[30] = {
		30, 0, 1, 29, 2, 9, 25, 10, 24, 12, 20, 13, 19, 15, 16, 31, 36};
	knownDFill[31] = {
		31, 2, 30, 3, 29, 9, 26, 11, 24, 12, 21, 14, 19, 15, 17, 32, 37};
	knownDFill[32] = {32, 2, 31, 4, 29, 9, 27, 12, 22, 15, 18, 33, 38};
	knownDFill[33] = {33, 2, 32, 5, 29, 9, 28, 12, 23, 15, 34, 39};
	knownDFill[34] = {34, 2, 33, 6, 29, 9, 12, 15, 40};
	knownDFill[35] = {
		35, 0, 36, 1, 3, 4, 5, 6, 7, 9, 29, 10, 11, 12, 24, 13, 14, 15, 19, 42};
	knownDFill[36] = {
		36, 0, 1, 35, 2, 9, 30, 10, 29, 12, 25, 13, 24, 15, 20, 16, 19, 37, 43};
	knownDFill[37] = {37, 2, 36, 3, 35, 9, 31, 11, 29, 12, 26, 14, 24, 15, 21,
		17, 19, 38, 44};
	knownDFill[38] = {38, 2, 37, 4, 35, 9, 32, 12, 27, 15, 22, 18, 19, 39, 45};
	knownDFill[39] = {39, 2, 38, 5, 35, 9, 33, 12, 28, 15, 23, 40, 46};
	knownDFill[40] = {40, 2, 39, 6, 35, 9, 34, 12, 41, 47};
	knownDFill[41] = {41, 2, 40, 7, 35, 9, 12, 48};
	knownDFill[42] = {42, 0, 43, 1, 9, 35, 12, 29, 15, 24, 19};
	knownDFill[43] = {43, 0, 1, 42, 2, 9, 36, 10, 35, 12, 30, 13, 29, 15, 25,
		16, 24, 19, 20, 44, 49};
	knownDFill[44] = {44, 2, 43, 9, 37, 11, 35, 12, 31, 14, 29, 15, 26, 17, 24,
		19, 21, 45, 50};
	knownDFill[45] = {45, 2, 44, 9, 38, 12, 32, 15, 27, 18, 24, 19, 22, 46, 51};
	knownDFill[46] = {46, 2, 45, 9, 39, 12, 33, 15, 28, 19, 23, 47, 52};
	knownDFill[47] = {47, 2, 46, 9, 40, 12, 34, 48, 53};
	knownDFill[48] = {48, 2, 47, 9, 41, 54};
	knownDFill[49] = {
		49, 2, 9, 42, 43, 12, 35, 36, 13, 15, 29, 30, 16, 19, 24, 25, 20, 50};
	knownDFill[50] = {
		50, 2, 49, 9, 44, 12, 37, 14, 35, 15, 31, 17, 29, 19, 26, 21, 24, 51};
	knownDFill[51] = {
		51, 2, 50, 9, 45, 12, 38, 15, 32, 18, 29, 19, 27, 22, 24, 52};
	knownDFill[52] = {52, 2, 51, 9, 46, 12, 39, 15, 33, 19, 28, 23, 24, 53};
	knownDFill[53] = {53, 2, 52, 9, 47, 12, 40, 15, 34, 54};
	knownDFill[54] = {54, 2, 53, 9, 48, 12, 41, 55};
	knownDFill[55] = {55, 2, 54};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 863);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 290600199485, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::He, 2), 44.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::V, 1), 46.0, 0.01);
	BOOST_REQUIRE_CLOSE(network.getTotalRadiusConcentration(dConcs, Spec::I, 1),
		0.15785478, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::T, 1), 99.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 2905985128900, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(
			fluxes[i], xolotl::test::TritiumSpeciesList::knownFluxes[i], 0.01);
	}

	// Check the partials computation
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	int startingIdx = 0;
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
				auto index = std::distance(row.begin(), iter);
				BOOST_REQUIRE_CLOSE(hPartials[startingIdx + index],
					xolotl::test::TritiumSpeciesList::knownPartials
						[startingIdx + j],
					0.01);
			}
			startingIdx += row.size();
		}
	}

	// Check clusters
	NetworkType::Composition comp = NetworkType::Composition::zero();
	comp[Spec::V] = 1;
	auto cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 1);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.137265, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 3.6, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 505312.69, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 1800000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 1.3, 0.01);
	const auto& reg0 = cluster.getRegion();
	Composition lo = reg0.getOrigin();
	Composition hi = reg0.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 1);
	auto momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);

	comp[Spec::V] = 0;
	comp[Spec::I] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 0);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.15785, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 10.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 78358278338, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 88000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.01, 0.01);
	const auto& reg1 = cluster.getRegion();
	lo = reg1.getOrigin();
	hi = reg1.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);

	comp[Spec::I] = 0;
	comp[Spec::He] = 5;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 24);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.3648, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getFormationEnergy(), 26.1, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 1242214406, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 5000000000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.12, 0.01);
	const auto& reg2 = cluster.getRegion();
	lo = reg2.getOrigin();
	hi = reg2.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 5);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 6);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 1);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);

	comp[Spec::He] = 0;
	comp[Spec::T] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 2);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.075, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getFormationEnergy(), std::numeric_limits<double>::infinity());
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		cluster.getDiffusionCoefficient(gridId), 2808640007, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 231000000000, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getMigrationEnergy(), 0.38, 0.01);
	const auto& reg3 = cluster.getRegion();
	lo = reg3.getOrigin();
	hi = reg3.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 2);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);

	comp[Spec::He] = 3;
	comp[Spec::V] = 1;
	cluster = network.findCluster(comp, plsm::HostMemSpace{});
	BOOST_REQUIRE_EQUAL(cluster.getId(), 17);
	BOOST_REQUIRE_CLOSE(cluster.getReactionRadius(), 0.137265, 0.01);
	BOOST_REQUIRE_EQUAL(cluster.getFormationEnergy(), 11.5304);
	BOOST_REQUIRE_CLOSE(cluster.getTemperature(gridId), 1000.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionCoefficient(gridId), 0.0, 0.01);
	BOOST_REQUIRE_CLOSE(cluster.getDiffusionFactor(), 0.0, 0.01);
	BOOST_REQUIRE_EQUAL(
		cluster.getMigrationEnergy(), std::numeric_limits<double>::infinity());
	const auto& reg4 = cluster.getRegion();
	lo = reg4.getOrigin();
	hi = reg4.getUpperLimitPoint();
	BOOST_REQUIRE_EQUAL(lo[Spec::V], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::V], 2);
	BOOST_REQUIRE_EQUAL(lo[Spec::I], 0);
	BOOST_REQUIRE_EQUAL(hi[Spec::I], 1);
	BOOST_REQUIRE_EQUAL(lo[Spec::He], 3);
	BOOST_REQUIRE_EQUAL(hi[Spec::He], 4);
	BOOST_REQUIRE_EQUAL(lo[Spec::T], 1);
	BOOST_REQUIRE_EQUAL(hi[Spec::T], 2);
	momId = cluster.getMomentIds();
	BOOST_REQUIRE_EQUAL(momId.extent(0), 3);
}

BOOST_AUTO_TEST_CASE(HeVGrouped)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 0 50 6" << std::endl
			  << "process=reaction" << std::endl
			  << "material=W100" << std::endl
			  << "grouping=31 4 4" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIHeliumSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	auto network = dynamic_pointer_cast<NetworkType>(
		factory::network::NetworkHandlerFactory::get()
			.generate(opts)
			->getNetwork());

	BOOST_REQUIRE_EQUAL(network->getNumClusters(), 2874);
	BOOST_REQUIRE_EQUAL(network->getDOF(), 3300);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network->getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 2874);
	typename NetworkType::PhaseSpace phaseSpace = network->getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 3);

	BOOST_REQUIRE_EQUAL(network->getNumberOfSpecies(), 3);
	BOOST_REQUIRE_EQUAL(network->getNumberOfSpeciesNoI(), 2);
}

BOOST_AUTO_TEST_CASE(HeDVGrouped)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 1 0 20 6" << std::endl
			  << "process=reaction" << std::endl
			  << "material=W100" << std::endl
			  << "grouping=10 4 4" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIDeuteriumSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	auto network = dynamic_pointer_cast<NetworkType>(
		factory::network::NetworkHandlerFactory::get()
			.generate(opts)
			->getNetwork());

	BOOST_REQUIRE_EQUAL(network->getNumClusters(), 4967);
	BOOST_REQUIRE_EQUAL(network->getDOF(), 6656);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network->getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 4967);
	typename NetworkType::PhaseSpace phaseSpace = network->getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 4);

	BOOST_REQUIRE_EQUAL(network->getNumberOfSpecies(), 4);
	BOOST_REQUIRE_EQUAL(network->getNumberOfSpeciesNoI(), 3);
}

BOOST_AUTO_TEST_CASE(HeTVGrouped)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 1 20 6" << std::endl
			  << "process=reaction" << std::endl
			  << "material=W100" << std::endl
			  << "grouping=10 4 4" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSITritiumSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	auto network = dynamic_pointer_cast<NetworkType>(
		factory::network::NetworkHandlerFactory::get()
			.generate(opts)
			->getNetwork());

	BOOST_REQUIRE_EQUAL(network->getNumClusters(), 4967);
	BOOST_REQUIRE_EQUAL(network->getDOF(), 6656);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network->getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 4967);
	typename NetworkType::PhaseSpace phaseSpace = network->getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 4);

	BOOST_REQUIRE_EQUAL(network->getNumberOfSpecies(), 4);
	BOOST_REQUIRE_EQUAL(network->getNumberOfSpeciesNoI(), 3);
}

BOOST_AUTO_TEST_CASE(HeDTVGrouped)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 1 1 6 2" << std::endl
			  << "process=reaction" << std::endl
			  << "material=W100" << std::endl
			  << "grouping=3 2 2" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIFullSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	auto network = dynamic_pointer_cast<NetworkType>(
		factory::network::NetworkHandlerFactory::get()
			.generate(opts)
			->getNetwork());

	BOOST_REQUIRE_EQUAL(network->getNumClusters(), 2383);
	BOOST_REQUIRE_EQUAL(network->getDOF(), 4259);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network->getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 2383);
	typename NetworkType::PhaseSpace phaseSpace = network->getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 5);

	BOOST_REQUIRE_EQUAL(network->getNumberOfSpecies(), 5);
	BOOST_REQUIRE_EQUAL(network->getNumberOfSpeciesNoI(), 4);
}

BOOST_AUTO_TEST_CASE(IGrouped)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=0 0 0 100 10000" << std::endl
			  << "process=reaction" << std::endl
			  << "material=W100" << std::endl
			  << "grouping=101 2 2" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIHeliumSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	auto network = dynamic_pointer_cast<NetworkType>(
		factory::network::NetworkHandlerFactory::get()
			.generate(opts)
			->getNetwork());

	BOOST_REQUIRE_EQUAL(network->getNumClusters(), 813);
	BOOST_REQUIRE_EQUAL(network->getDOF(), 1424);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network->getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 813);
	typename NetworkType::PhaseSpace phaseSpace = network->getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 3);

	BOOST_REQUIRE_EQUAL(network->getNumberOfSpecies(), 3);
	BOOST_REQUIRE_EQUAL(network->getNumberOfSpeciesNoI(), 2);
}

BOOST_AUTO_TEST_CASE(LargeBubble)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=8 0 0 5 2" << std::endl
			  << "process=reaction largeBubble" << std::endl
			  << "tempParam=1000" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIHeliumSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	// Get the temperature and lattice constant from the options
	double latticeConst = opts.getLatticeParameter() <= 0.0 ?
		xolotl::core::tungstenLatticeConstant :
		opts.getLatticeParameter();
	double temperature = opts.getTempParam();

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe =
		util::getMaxHePerVLoop(maxV, latticeConst, temperature);
	NetworkType network({maxHe, maxV, maxI}, 1, opts);

	BOOST_REQUIRE(!network.hasDeuterium());
	BOOST_REQUIRE(!network.hasTritium());

	BOOST_REQUIRE_EQUAL(network.getNumClusters(), 103);
	BOOST_REQUIRE_EQUAL(network.getDOF(), 106);

	BOOST_REQUIRE_CLOSE(network.getLatticeParameter(), 0.317, 0.01);
	BOOST_REQUIRE_CLOSE(network.getAtomicVolume(), 0.0159275, 0.01);
	BOOST_REQUIRE_CLOSE(network.getImpurityRadius(), 0.3, 0.01);
	BOOST_REQUIRE_CLOSE(network.getInterstitialBias(), 1.15, 0.01);

	BOOST_REQUIRE(network.getEnableStdReaction() == true);
	BOOST_REQUIRE(network.getEnableLargeBubble() == true);

	BOOST_REQUIRE_EQUAL(network.getGridSize(), 1);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network.getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 103);
	typename NetworkType::PhaseSpace phaseSpace = network.getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 3);

	BOOST_REQUIRE_EQUAL(network.getNumberOfSpecies(), 3);
	BOOST_REQUIRE_EQUAL(network.getNumberOfSpeciesNoI(), 2);

	// Get the diagonal fill
	const auto dof = network.getDOF();
	NetworkType::SparseFillMap knownDFill;
	knownDFill[0] = {0, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 20,
		21, 22, 23, 24, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 38, 39, 40, 41,
		42, 44, 45, 46, 47, 48, 50, 51, 52, 53, 54, 56, 57, 58, 59, 61, 62, 63,
		65, 66, 67, 69, 70, 71, 73, 74, 75, 77, 78, 79, 81, 82, 84, 85, 87, 88,
		90, 91, 93, 95, 1, 7, 103, 104, 105, 49, 55, 76, 89, 94, 102, 13, 43,
		72, 86, 92, 101, 19, 37, 68, 83, 100, 25, 31, 64, 80, 99, 60, 98, 97,
		96};
	knownDFill[1] = {1, 0, 2, 3, 4, 5, 6, 9, 10, 11, 12, 15, 16, 17, 18, 21, 22,
		23, 24, 27, 28, 29, 30, 33, 34, 35, 36, 39, 40, 41, 42, 45, 46, 47, 48,
		51, 52, 53, 54, 57, 58, 59, 62, 63, 66, 67, 70, 71, 74, 75, 78, 79, 82,
		85, 88, 91, 13, 49, 19, 43, 89, 25, 37, 86, 102, 31, 76, 83, 101, 55,
		72, 80, 100, 50, 68, 77, 99};
	knownDFill[2] = {2, 0, 3, 1, 4, 103, 104, 105, 5, 6, 7, 8, 9, 10, 11, 12,
		13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
		31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
		49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
		67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
		85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
		102};
	knownDFill[3] = {3, 0, 4, 1, 5, 2, 7, 13, 19, 25, 31, 37, 43, 9};
	knownDFill[4] = {4, 0, 5, 1, 6, 2, 3, 7, 13, 19, 25, 31, 37, 43, 10};
	knownDFill[5] = {5, 0, 6, 1, 2, 4, 7, 13, 19, 25, 31, 37, 43, 11};
	knownDFill[6] = {6, 0, 1, 2, 5, 7, 13, 19, 25, 31, 37, 43, 12};
	knownDFill[7] = {7, 0, 8, 1, 9, 2, 3, 4, 5, 6, 103, 104, 105, 10, 11, 12,
		13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
		31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
		49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
		67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
		85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
		102};
	knownDFill[8] = {8, 0, 9, 1, 10, 2, 7, 13, 19, 25, 31, 37, 43, 14};
	knownDFill[9] = {9, 0, 10, 1, 11, 2, 8, 3, 7, 13, 19, 25, 31, 37, 43, 15};
	knownDFill[10] = {10, 0, 11, 1, 12, 2, 9, 4, 7, 13, 19, 25, 31, 37, 43, 16};
	knownDFill[11] = {11, 0, 12, 1, 2, 10, 5, 7, 13, 19, 25, 31, 37, 43, 17};
	knownDFill[12] = {12, 0, 1, 2, 11, 6, 7, 13, 19, 25, 31, 37, 43, 18};
	knownDFill[13] = {13, 0, 14, 1, 15, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 103,
		104, 105, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
		31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
		49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
		67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
		85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
		102};
	knownDFill[14] = {14, 0, 15, 1, 16, 2, 13, 7, 8, 19, 25, 31, 37, 43, 20};
	knownDFill[15] = {
		15, 0, 16, 1, 17, 2, 14, 3, 13, 7, 9, 19, 25, 31, 37, 43, 21};
	knownDFill[16] = {
		16, 0, 17, 1, 18, 2, 15, 4, 13, 7, 10, 19, 25, 31, 37, 43, 22};
	knownDFill[17] = {
		17, 0, 18, 1, 2, 16, 5, 13, 7, 11, 19, 25, 31, 37, 43, 23};
	knownDFill[18] = {18, 0, 1, 2, 17, 6, 13, 7, 12, 19, 25, 31, 37, 43, 24};
	knownDFill[19] = {19, 0, 20, 1, 21, 2, 3, 4, 5, 6, 7, 13, 8, 9, 10, 11, 12,
		14, 15, 16, 17, 18, 103, 104, 105, 22, 23, 24, 25, 26, 27, 28, 29, 30,
		31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
		49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
		67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
		85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
		102};
	knownDFill[20] = {
		20, 0, 21, 1, 22, 2, 19, 7, 14, 8, 13, 25, 31, 37, 43, 26};
	knownDFill[21] = {
		21, 0, 22, 1, 23, 2, 20, 3, 19, 7, 15, 9, 13, 25, 31, 37, 43, 27};
	knownDFill[22] = {
		22, 0, 23, 1, 24, 2, 21, 4, 19, 7, 16, 10, 13, 25, 31, 37, 43, 28};
	knownDFill[23] = {
		23, 0, 24, 1, 2, 22, 5, 19, 7, 17, 11, 13, 25, 31, 37, 43, 29};
	knownDFill[24] = {
		24, 0, 1, 2, 23, 6, 19, 7, 18, 12, 13, 25, 31, 37, 43, 30};
	knownDFill[25] = {25, 0, 26, 1, 27, 2, 3, 4, 5, 6, 7, 19, 8, 9, 10, 11, 12,
		13, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 103, 104, 105, 28, 29, 30,
		31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
		49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
		67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
		85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
		102};
	knownDFill[26] = {
		26, 0, 27, 1, 28, 2, 25, 7, 20, 8, 19, 13, 14, 31, 37, 43, 32};
	knownDFill[27] = {
		27, 0, 28, 1, 29, 2, 26, 3, 25, 7, 21, 9, 19, 13, 15, 31, 37, 43, 33};
	knownDFill[28] = {
		28, 0, 29, 1, 30, 2, 27, 4, 25, 7, 22, 10, 19, 13, 16, 31, 37, 43, 34};
	knownDFill[29] = {
		29, 0, 30, 1, 2, 28, 5, 25, 7, 23, 11, 19, 13, 17, 31, 37, 43, 35};
	knownDFill[30] = {
		30, 0, 1, 2, 29, 6, 25, 7, 24, 12, 19, 13, 18, 31, 37, 43, 36};
	knownDFill[31] = {31, 0, 32, 1, 33, 2, 3, 4, 5, 6, 7, 25, 8, 9, 10, 11, 12,
		13, 19, 14, 15, 16, 17, 18, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 103,
		104, 105, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48,
		49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
		67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
		85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
		102};
	knownDFill[32] = {
		32, 0, 33, 1, 34, 2, 31, 7, 26, 8, 25, 13, 20, 14, 19, 37, 43, 38};
	knownDFill[33] = {33, 0, 34, 1, 35, 2, 32, 3, 31, 7, 27, 9, 25, 13, 21, 15,
		19, 37, 43, 39};
	knownDFill[34] = {34, 0, 35, 1, 36, 2, 33, 4, 31, 7, 28, 10, 25, 13, 22, 16,
		19, 37, 43, 40};
	knownDFill[35] = {
		35, 0, 36, 1, 2, 34, 5, 31, 7, 29, 11, 25, 13, 23, 17, 19, 37, 43, 41};
	knownDFill[36] = {
		36, 0, 1, 2, 35, 6, 31, 7, 30, 12, 25, 13, 24, 18, 19, 37, 43, 42};
	knownDFill[37] = {37, 0, 38, 1, 39, 2, 3, 4, 5, 6, 7, 31, 8, 9, 10, 11, 12,
		13, 25, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30,
		32, 33, 34, 35, 36, 103, 104, 105, 40, 41, 42, 43, 44, 45, 46, 47, 48,
		49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66,
		67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84,
		85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101,
		102};
	knownDFill[38] = {
		38, 0, 39, 1, 40, 2, 37, 7, 32, 8, 31, 13, 26, 14, 25, 19, 20, 43, 44};
	knownDFill[39] = {39, 0, 40, 1, 41, 2, 38, 3, 37, 7, 33, 9, 31, 13, 27, 15,
		25, 19, 21, 43, 45};
	knownDFill[40] = {40, 0, 41, 1, 42, 2, 39, 4, 37, 7, 34, 10, 31, 13, 28, 16,
		25, 19, 22, 43, 46};
	knownDFill[41] = {41, 0, 42, 1, 2, 40, 5, 37, 7, 35, 11, 31, 13, 29, 17, 25,
		19, 23, 43, 47};
	knownDFill[42] = {
		42, 0, 1, 2, 41, 6, 37, 7, 36, 12, 31, 13, 30, 18, 25, 19, 24, 43, 48};
	knownDFill[43] = {43, 0, 44, 1, 45, 2, 3, 4, 5, 6, 7, 37, 8, 9, 10, 11, 12,
		13, 31, 14, 15, 16, 17, 18, 19, 25, 20, 21, 22, 23, 24, 26, 27, 28, 29,
		30, 32, 33, 34, 35, 36, 38, 39, 40, 41, 42, 103, 104, 105, 46, 47, 48,
		50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67,
		68, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86,
		87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 49};
	knownDFill[44] = {44, 0, 45, 1, 46, 2, 43, 7, 38, 8, 37, 13, 32, 14, 31, 19,
		26, 20, 25, 50};
	knownDFill[45] = {45, 0, 46, 1, 47, 2, 44, 3, 43, 7, 39, 9, 37, 13, 33, 15,
		31, 19, 27, 21, 25, 51};
	knownDFill[46] = {46, 0, 47, 1, 48, 2, 45, 4, 43, 7, 40, 10, 37, 13, 34, 16,
		31, 19, 28, 22, 25, 52};
	knownDFill[47] = {47, 0, 48, 1, 2, 46, 5, 43, 7, 41, 11, 37, 13, 35, 17, 31,
		19, 29, 23, 25, 53};
	knownDFill[48] = {48, 0, 1, 2, 47, 6, 43, 7, 42, 12, 37, 13, 36, 18, 31, 19,
		30, 24, 25, 54};
	knownDFill[49] = {49, 0, 50, 1, 51, 2, 7, 43, 13, 37, 19, 31, 25};
	knownDFill[50] = {50, 0, 51, 1, 52, 2, 49, 7, 44, 8, 43, 13, 38, 14, 37, 19,
		32, 20, 31, 25, 26, 55};
	knownDFill[51] = {51, 0, 52, 1, 53, 2, 50, 7, 45, 9, 43, 13, 39, 15, 37, 19,
		33, 21, 31, 25, 27, 56};
	knownDFill[52] = {52, 0, 53, 1, 54, 2, 51, 7, 46, 10, 43, 13, 40, 16, 37,
		19, 34, 22, 31, 25, 28, 57};
	knownDFill[53] = {53, 0, 54, 1, 2, 52, 7, 47, 11, 43, 13, 41, 17, 37, 19,
		35, 23, 31, 25, 29, 58};
	knownDFill[54] = {54, 0, 1, 2, 53, 7, 48, 12, 43, 13, 42, 18, 37, 19, 36,
		24, 31, 25, 30, 59};
	knownDFill[55] = {55, 0, 56, 1, 57, 2, 7, 49, 50, 13, 43, 44, 14, 19, 37,
		38, 20, 25, 31, 32, 26};
	knownDFill[56] = {56, 0, 57, 1, 58, 2, 55, 7, 51, 13, 45, 15, 43, 19, 39,
		21, 37, 25, 33, 27, 31, 60};
	knownDFill[57] = {57, 0, 58, 1, 59, 2, 56, 7, 52, 13, 46, 16, 43, 19, 40,
		22, 37, 25, 34, 28, 31, 61};
	knownDFill[58] = {58, 0, 59, 1, 2, 57, 7, 53, 13, 47, 17, 43, 19, 41, 23,
		37, 25, 35, 29, 31, 62};
	knownDFill[59] = {59, 0, 1, 2, 58, 7, 54, 13, 48, 18, 43, 19, 42, 24, 37,
		25, 36, 30, 31, 63};
	knownDFill[60] = {60, 0, 61, 1, 62, 2, 7, 55, 56, 13, 49, 50, 51, 19, 43,
		44, 45, 20, 21, 25, 37, 38, 39, 26, 27, 31, 32, 33, 64};
	knownDFill[61] = {61, 0, 62, 1, 63, 2, 60, 7, 57, 13, 52, 19, 46, 22, 43,
		25, 40, 28, 37, 31, 34, 65};
	knownDFill[62] = {62, 0, 63, 1, 2, 61, 7, 58, 13, 53, 19, 47, 23, 43, 25,
		41, 29, 37, 31, 35, 66};
	knownDFill[63] = {63, 0, 1, 2, 62, 7, 59, 13, 54, 19, 48, 24, 43, 25, 42,
		30, 37, 31, 36, 67};
	knownDFill[64] = {64, 0, 65, 1, 66, 2, 7, 60, 13, 55, 56, 19, 49, 50, 51,
		25, 43, 44, 45, 26, 27, 31, 37, 38, 39, 32, 33, 68};
	knownDFill[65] = {65, 0, 66, 1, 67, 2, 64, 7, 61, 13, 57, 19, 52, 25, 46,
		28, 43, 31, 40, 34, 37, 69};
	knownDFill[66] = {66, 0, 67, 1, 2, 65, 7, 62, 13, 58, 19, 53, 25, 47, 29,
		43, 31, 41, 35, 37, 70};
	knownDFill[67] = {67, 0, 1, 2, 66, 7, 63, 13, 59, 19, 54, 25, 48, 30, 43,
		31, 42, 36, 37, 71};
	knownDFill[68] = {68, 0, 69, 1, 70, 2, 7, 64, 13, 60, 19, 55, 56, 25, 49,
		50, 51, 31, 43, 44, 45, 32, 33, 37, 38, 39, 72};
	knownDFill[69] = {69, 0, 70, 1, 71, 2, 68, 7, 65, 13, 61, 19, 57, 25, 52,
		31, 46, 34, 43, 37, 40, 73};
	knownDFill[70] = {70, 0, 71, 1, 2, 69, 7, 66, 13, 62, 19, 58, 25, 53, 31,
		47, 35, 43, 37, 41, 74};
	knownDFill[71] = {71, 0, 1, 2, 70, 7, 67, 13, 63, 19, 59, 25, 54, 31, 48,
		36, 43, 37, 42, 75};
	knownDFill[72] = {72, 0, 73, 1, 74, 2, 7, 68, 13, 64, 19, 60, 25, 55, 56,
		31, 49, 50, 51, 37, 43, 44, 45, 38, 39, 76};
	knownDFill[73] = {73, 0, 74, 1, 75, 2, 72, 7, 69, 13, 65, 19, 61, 25, 57,
		31, 52, 37, 46, 40, 43, 77};
	knownDFill[74] = {74, 0, 75, 1, 2, 73, 7, 70, 13, 66, 19, 62, 25, 58, 31,
		53, 37, 47, 41, 43, 78};
	knownDFill[75] = {75, 0, 1, 2, 74, 7, 71, 13, 67, 19, 63, 25, 59, 31, 54,
		37, 48, 42, 43, 79};
	knownDFill[76] = {76, 0, 77, 1, 78, 2, 7, 72, 13, 68, 19, 64, 25, 60, 31,
		55, 56, 37, 49, 50, 51, 43, 44, 45};
	knownDFill[77] = {77, 0, 78, 1, 79, 2, 76, 7, 73, 13, 69, 19, 65, 25, 61,
		31, 57, 37, 52, 43, 46, 80};
	knownDFill[78] = {78, 0, 79, 1, 2, 77, 7, 74, 13, 70, 19, 66, 25, 62, 31,
		58, 37, 53, 43, 47, 81};
	knownDFill[79] = {79, 0, 1, 2, 78, 7, 75, 13, 71, 19, 67, 25, 63, 31, 59,
		37, 54, 43, 48, 82};
	knownDFill[80] = {80, 0, 81, 1, 82, 2, 7, 76, 77, 13, 72, 73, 19, 68, 69,
		25, 64, 65, 31, 60, 61, 37, 55, 56, 57, 43, 50, 51, 52, 83};
	knownDFill[81] = {81, 0, 82, 2, 80, 7, 78, 13, 74, 19, 70, 25, 66, 31, 62,
		37, 58, 43, 53, 84};
	knownDFill[82] = {82, 0, 1, 2, 81, 7, 79, 13, 75, 19, 71, 25, 67, 31, 63,
		37, 59, 43, 54, 85};
	knownDFill[83] = {83, 0, 84, 1, 85, 2, 7, 80, 13, 76, 77, 19, 72, 73, 25,
		68, 69, 31, 64, 65, 37, 60, 61, 43, 55, 56, 57, 86};
	knownDFill[84] = {84, 0, 85, 2, 83, 7, 81, 13, 78, 19, 74, 25, 70, 31, 66,
		37, 62, 43, 58, 87};
	knownDFill[85] = {85, 0, 1, 2, 84, 7, 82, 13, 79, 19, 75, 25, 71, 31, 67,
		37, 63, 43, 59, 88};
	knownDFill[86] = {86, 0, 87, 1, 88, 2, 7, 83, 13, 80, 19, 76, 77, 25, 72,
		73, 31, 68, 69, 37, 64, 65, 43, 60, 61, 89};
	knownDFill[87] = {87, 0, 88, 2, 86, 7, 84, 13, 81, 19, 78, 25, 74, 31, 70,
		37, 66, 43, 62, 90};
	knownDFill[88] = {88, 0, 1, 2, 87, 7, 85, 13, 82, 19, 79, 25, 75, 31, 71,
		37, 67, 43, 63, 91};
	knownDFill[89] = {89, 0, 90, 1, 91, 2, 7, 86, 13, 83, 19, 80, 25, 76, 77,
		31, 72, 73, 37, 68, 69, 43, 64, 65};
	knownDFill[90] = {90, 0, 91, 2, 89, 7, 87, 13, 84, 19, 81, 25, 78, 31, 74,
		37, 70, 43, 66, 92};
	knownDFill[91] = {91, 0, 1, 2, 90, 7, 88, 13, 85, 19, 82, 25, 79, 31, 75,
		37, 71, 43, 67, 93};
	knownDFill[92] = {92, 0, 93, 2, 7, 89, 90, 13, 86, 87, 19, 83, 84, 25, 80,
		81, 31, 76, 77, 78, 37, 72, 73, 74, 43, 68, 69, 70, 94};
	knownDFill[93] = {93, 0, 2, 92, 7, 91, 13, 88, 19, 85, 25, 82, 31, 79, 37,
		75, 43, 71, 95};
	knownDFill[94] = {94, 0, 95, 2, 7, 92, 13, 89, 90, 19, 86, 87, 25, 83, 84,
		31, 80, 81, 37, 76, 77, 78, 43, 72, 73, 74};
	knownDFill[95] = {95, 0, 2, 94, 7, 93, 13, 91, 19, 88, 25, 85, 31, 82, 37,
		79, 43, 75, 96};
	knownDFill[96] = {96, 2, 7, 94, 95, 13, 92, 93, 19, 89, 90, 91, 25, 86, 87,
		88, 31, 83, 84, 85, 37, 80, 81, 82, 43, 77, 78, 79, 97};
	knownDFill[97] = {97, 2, 7, 96, 13, 94, 95, 19, 92, 93, 25, 89, 90, 91, 31,
		86, 87, 88, 37, 83, 84, 85, 43, 80, 81, 82, 98};
	knownDFill[98] = {98, 2, 7, 97, 13, 96, 19, 94, 95, 25, 92, 93, 31, 89, 90,
		91, 37, 86, 87, 88, 43, 83, 84, 85, 99};
	knownDFill[99] = {99, 2, 7, 98, 13, 97, 19, 96, 25, 94, 95, 31, 92, 93, 37,
		89, 90, 91, 43, 86, 87, 88, 100};
	knownDFill[100] = {100, 2, 7, 99, 13, 98, 19, 97, 25, 96, 31, 94, 95, 37,
		92, 93, 43, 89, 90, 91, 101};
	knownDFill[101] = {101, 2, 7, 100, 13, 99, 19, 98, 25, 97, 31, 96, 37, 94,
		95, 43, 92, 93, 102};
	knownDFill[102] = {
		102, 2, 7, 101, 13, 100, 19, 99, 25, 98, 31, 97, 37, 96, 43, 94, 95};
	knownDFill[103] = {103, 2, 104, 105, 6, 12, 18, 24, 30, 36, 42, 48, 54, 59,
		63, 67, 71, 75, 79, 82, 85, 88, 91, 93, 95, 96, 97, 98, 99, 100, 101,
		102, 7, 13, 19, 25, 31, 37, 43};
	knownDFill[104] = {104, 2, 103, 105, 6, 12, 18, 24, 30, 36, 42, 48, 54, 59,
		63, 67, 71, 75, 79, 82, 85, 88, 91, 93, 95, 96, 97, 98, 99, 100, 101,
		102, 7, 13, 19, 25, 31, 37, 43};
	knownDFill[105] = {105, 2, 103, 104, 6, 12, 18, 24, 30, 36, 42, 48, 54, 59,
		63, 67, 71, 75, 79, 82, 85, 88, 91, 93, 95, 96, 97, 98, 99, 100, 101,
		102, 7, 13, 19, 25, 31, 37, 43};

	NetworkType::SparseFillMap dfill;
	auto nPartials = network.getDiagonalFill(dfill);
	BOOST_REQUIRE_EQUAL(nPartials, 3030);
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			BOOST_REQUIRE_EQUAL(row.size(), knownDFill[i].size());
		}
	}

	// Set temperatures
	std::vector<double> temperatures = {1000.0};
	std::vector<double> depths = {1.0};
	network.setTemperatures(temperatures, depths);
	NetworkType::IndexType gridId = 0;

	// Check the largest rate
	BOOST_REQUIRE_CLOSE(network.getLargestRate(), 621745254537, 0.01);

	// Create a concentration vector where every field is at 1.0
	std::vector<double> concentrations(dof + 1, 1.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hConcs = HostUnmanaged(concentrations.data(), dof + 1);
	auto dConcs = Kokkos::View<double*>("Concentrations", dof + 1);
	deep_copy(dConcs, hConcs);

	// Check the total concentrations
	BOOST_REQUIRE_CLOSE(
		network.getTotalConcentration(dConcs, Spec::He, 2), 90.0, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalAtomConcentration(dConcs, Spec::V, 1), 321.0, 0.01);
	BOOST_REQUIRE_CLOSE(network.getTotalRadiusConcentration(dConcs, Spec::I, 1),
		0.356278, 0.01);
	BOOST_REQUIRE_CLOSE(
		network.getTotalTrappedAtomConcentration(dConcs, Spec::He, 1), 909.0,
		0.01);

	// Check the left side rate of 0th cluster
	BOOST_REQUIRE_CLOSE(
		network.getLeftSideRate(dConcs, 0, gridId), 27551779352980, 0.01);

	// Create a flux vector where every field is at 0.0
	std::vector<double> fluxes(dof + 1, 0.0);
	using HostUnmanaged =
		Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
	auto hFluxes = HostUnmanaged(fluxes.data(), dof + 1);
	auto dFluxes = Kokkos::View<double*>("Fluxes", dof + 1);
	deep_copy(dFluxes, hFluxes);

	// Check the fluxes computation
	network.computeAllFluxes(dConcs, dFluxes, gridId);
	deep_copy(hFluxes, dFluxes);
	for (NetworkType::IndexType i = 0; i < dof + 1; i++) {
		BOOST_REQUIRE_CLOSE(
			fluxes[i], xolotl::test::LargeBubble::knownFluxes[i], 0.01);
	}

	// Check the partials computation
	auto vals = Kokkos::View<double*>("solverPartials", nPartials);
	network.computeAllPartials(dConcs, vals, gridId);
	auto hPartials = create_mirror_view(vals);
	deep_copy(hPartials, vals);
	int startingIdx = 0;
	for (NetworkType::IndexType i = 0; i < dof; i++) {
		auto rowIter = dfill.find(i);
		if (rowIter != dfill.end()) {
			const auto& row = rowIter->second;
			for (NetworkType::IndexType j = 0; j < row.size(); j++) {
				auto iter = find(row.begin(), row.end(), knownDFill[i][j]);
				auto index = std::distance(row.begin(), iter);
				BOOST_REQUIRE_CLOSE(hPartials[startingIdx + index],
					xolotl::test::LargeBubble::knownPartials[startingIdx + j],
					0.01);
			}
			startingIdx += row.size();
		}
	}
}

BOOST_AUTO_TEST_CASE(VIGrouped)
{
	// Create the option to create a network
	xolotl::options::ConfOptions opts;
	// Create a good parameter file
	std::string parameterFile = "param.txt";
	std::ofstream paramFile(parameterFile);
	paramFile << "netParam=0 0 0 1000 10000" << std::endl
			  << "process=reaction" << std::endl
			  << "material=W100" << std::endl
			  << "grouping=101 2 2" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
	opts.readParams(cl.argc, cl.argv);

	std::remove(parameterFile.c_str());

	using NetworkType = PSIReactionNetwork<PSIHeliumSpeciesList>;
	using Spec = NetworkType::Species;
	using Composition = NetworkType::Composition;

	auto network = dynamic_pointer_cast<NetworkType>(
		factory::network::NetworkHandlerFactory::get()
			.generate(opts)
			->getNetwork());

	BOOST_REQUIRE_EQUAL(network->getNumClusters(), 1091);
	BOOST_REQUIRE_EQUAL(network->getDOF(), 1978);

	// TODO: Test each value explicitly?
	typename NetworkType::Bounds bounds = network->getAllClusterBounds();
	BOOST_REQUIRE_EQUAL(bounds.size(), 1091);
	typename NetworkType::PhaseSpace phaseSpace = network->getPhaseSpace();
	BOOST_REQUIRE_EQUAL(phaseSpace.size(), 3);

	BOOST_REQUIRE_EQUAL(network->getNumberOfSpecies(), 3);
	BOOST_REQUIRE_EQUAL(network->getNumberOfSpeciesNoI(), 2);
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(traits_Full, *utf::label("traits"))
{
	using NetworkType = PSIReactionNetwork<PSIFullSpeciesList>;
	using Spec = NetworkType::Species;

	static_assert(ReactionNetworkTraits<NetworkType>::numSpecies == 5);
	BOOST_REQUIRE_EQUAL(numberOfSpecies<Spec>(), 5);
	BOOST_REQUIRE_EQUAL(numberOfInterstitialSpecies<Spec>(), 1);
	BOOST_REQUIRE_EQUAL(numberOfVacancySpecies<Spec>(), 1);

	// toLabelString
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::He), "He");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::D), "D");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::T), "T");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::V), "V");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::I), "I");

	// toNameString
	BOOST_REQUIRE_EQUAL(toNameString(Spec::He), "Helium");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::D), "Deuterium");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::T), "Tritium");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::V), "Vacancy");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::I), "Interstitial");

	// parseSpeciesId / getSpeciesLabel / getSpeciesName
	auto network = NetworkType();
	auto sid = network.parseSpeciesId("He");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::He);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "He");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Helium");
	sid = network.parseSpeciesId("D");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::D);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "D");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Deuterium");
	sid = network.parseSpeciesId("T");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::T);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "T");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Tritium");
	sid = network.parseSpeciesId("V");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::V);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "V");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Vacancy");
	sid = network.parseSpeciesId("I");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::I);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "I");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Interstitial");

	// Define the sequence
	using GroupingRange = SpeciesForGrouping<Spec, 5>;

	// I goes on V
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::He), 0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::D), 1);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::T), 2);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::V), 3);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::I), 3);
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(traits_HeD, *utf::label("traits"))
{
	using NetworkType = PSIReactionNetwork<PSIDeuteriumSpeciesList>;
	using Spec = NetworkType::Species;

	static_assert(ReactionNetworkTraits<NetworkType>::numSpecies == 4);
	BOOST_REQUIRE_EQUAL(numberOfSpecies<Spec>(), 4);
	BOOST_REQUIRE_EQUAL(numberOfInterstitialSpecies<Spec>(), 1);
	BOOST_REQUIRE_EQUAL(numberOfVacancySpecies<Spec>(), 1);

	// toLabelString
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::He), "He");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::D), "D");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::V), "V");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::I), "I");

	// toNameString
	BOOST_REQUIRE_EQUAL(toNameString(Spec::He), "Helium");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::D), "Deuterium");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::V), "Vacancy");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::I), "Interstitial");

	// parseSpeciesId / getSpeciesLabel / getSpeciesName
	auto network = NetworkType();
	auto sid = network.parseSpeciesId("He");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::He);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "He");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Helium");
	sid = network.parseSpeciesId("D");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::D);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "D");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Deuterium");
	sid = network.parseSpeciesId("V");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::V);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "V");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Vacancy");
	sid = network.parseSpeciesId("I");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::I);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "I");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Interstitial");

	// Define the sequence
	using GroupingRange = SpeciesForGrouping<Spec, 4>;

	// I goes on V
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::He), 0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::D), 1);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::V), 2);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::I), 2);
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(traits_HeT, *utf::label("traits"))
{
	using NetworkType = PSIReactionNetwork<PSITritiumSpeciesList>;
	using Spec = NetworkType::Species;

	static_assert(ReactionNetworkTraits<NetworkType>::numSpecies == 4);
	BOOST_REQUIRE_EQUAL(numberOfSpecies<Spec>(), 4);
	BOOST_REQUIRE_EQUAL(numberOfInterstitialSpecies<Spec>(), 1);
	BOOST_REQUIRE_EQUAL(numberOfVacancySpecies<Spec>(), 1);

	// toLabelString
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::He), "He");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::T), "T");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::V), "V");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::I), "I");

	// toNameString
	BOOST_REQUIRE_EQUAL(toNameString(Spec::He), "Helium");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::T), "Tritium");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::V), "Vacancy");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::I), "Interstitial");

	// parseSpeciesId / getSpeciesLabel / getSpeciesName
	auto network = NetworkType();
	auto sid = network.parseSpeciesId("He");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::He);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "He");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Helium");
	sid = network.parseSpeciesId("T");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::T);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "T");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Tritium");
	sid = network.parseSpeciesId("V");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::V);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "V");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Vacancy");
	sid = network.parseSpeciesId("I");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::I);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "I");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Interstitial");

	// Define the sequence
	using GroupingRange = SpeciesForGrouping<Spec, 4>;

	// I goes on V
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::He), 0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::T), 1);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::V), 2);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::I), 2);
}

BOOST_AUTO_TEST_CASE_WITH_DECOR(traits_He, *utf::label("traits"))
{
	using NetworkType = PSIReactionNetwork<PSIHeliumSpeciesList>;
	using Spec = NetworkType::Species;

	static_assert(ReactionNetworkTraits<NetworkType>::numSpecies == 3);
	BOOST_REQUIRE_EQUAL(numberOfSpecies<Spec>(), 3);
	BOOST_REQUIRE_EQUAL(numberOfInterstitialSpecies<Spec>(), 1);
	BOOST_REQUIRE_EQUAL(numberOfVacancySpecies<Spec>(), 1);

	// toLabelString
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::He), "He");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::V), "V");
	BOOST_REQUIRE_EQUAL(toLabelString(Spec::I), "I");

	// toNameString
	BOOST_REQUIRE_EQUAL(toNameString(Spec::He), "Helium");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::V), "Vacancy");
	BOOST_REQUIRE_EQUAL(toNameString(Spec::I), "Interstitial");

	// parseSpeciesId / getSpeciesLabel / getSpeciesName
	auto network = NetworkType();
	auto sid = network.parseSpeciesId("He");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::He);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "He");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Helium");
	sid = network.parseSpeciesId("V");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::V);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "V");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Vacancy");
	sid = network.parseSpeciesId("I");
	BOOST_REQUIRE(sid.cast<Spec>() == Spec::I);
	BOOST_REQUIRE_EQUAL(network.getSpeciesLabel(sid), "I");
	BOOST_REQUIRE_EQUAL(network.getSpeciesName(sid), "Interstitial");

	// Define the sequence
	using GroupingRange = SpeciesForGrouping<Spec, 3>;

	// I goes on V
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::He), 0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::V), 1);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(Spec::I), 1);
}

BOOST_AUTO_TEST_SUITE_END()
