#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/unit_test.hpp>
#include <unordered_map>
#include <experimental/PSIReactionNetwork.h>

using namespace std;
using namespace xolotlCore;

class KokkosContext {
public:
	KokkosContext() {
		::Kokkos::initialize();
	}

	~KokkosContext() {
		::Kokkos::finalize();
	}
};
BOOST_GLOBAL_FIXTURE(KokkosContext);

/**
 * This suite is responsible for testing the experimental networks.
 */
BOOST_AUTO_TEST_SUITE(ExpNetwork_testSuite)

/**
 * This operation checks the creation of a PSI network
 */
BOOST_AUTO_TEST_CASE(PSINetwork) {

	using NetworkType =
	experimental::PSIReactionNetwork<experimental::PSIFullSpeciesList>;
	using Spec = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Create the option to create a network
	xolotlCore::Options opts;
	// Create a good parameter file
    //NetworkType rNetwork({31, 31, 31, 31, 31}, {{{2, 2, 2, 2, 2}}}, 0, Options{});
	std::ofstream paramFile("param.txt");
	paramFile << "netParam=8 1 1 6 6" << std::endl;
	paramFile.close();

	// Create a fake command line to read the options
	int argc = 2;
	char **argv = new char*[3];
	std::string appName = "fakeXolotlAppNameForTests";
	argv[0] = new char[appName.length() + 1];
	strcpy(argv[0], appName.c_str());
	std::string parameterFile = "param.txt";
	argv[1] = new char[parameterFile.length() + 1];
	strcpy(argv[1], parameterFile.c_str());
	argv[2] = 0; // null-terminate the array
	opts.readParams(argc, argv);

	// Get the boundaries from the options
	NetworkType::AmountType maxV = opts.getMaxV();
	NetworkType::AmountType maxI = opts.getMaxI();
	NetworkType::AmountType maxHe = experimental::PSIClusterGenerator<
			experimental::PSIFullSpeciesList>::getMaxHePerV(maxV);
	NetworkType::AmountType maxD = 2.0 / 3.0 * (double) maxHe;
	NetworkType::AmountType maxT = 2.0 / 3.0 * (double) maxHe;
	if (opts.getMaxImpurity() <= 0)
		maxHe = 0;
	if (opts.getMaxD() <= 0)
		maxD = 0;
	if (opts.getMaxT() <= 0)
		maxT = 0;
	NetworkType::SubdivisionRatio ratio = {31, 3, 3, 7, 7};
	std::vector<NetworkType::SubdivisionRatio> input;
	input.push_back(ratio);
	NetworkType::SubdivisionRatio ratioBis = {1, 7, 7, 1, 1};
	input.push_back(ratioBis);
	std::cout << maxHe << " " << maxD << " " << maxT << " " << maxV << " " << maxI << std::endl;
	NetworkType rNetwork( { maxHe, maxD, maxT, maxV, maxI }, input, 0, opts);

    rNetwork.syncClusterDataOnHost();
    rNetwork.getSubpaving().syncZones(plsm::onHost);

//	// PRINT ALL THE CLUSTERS
//	std::size_t numClusters = rNetwork.getDOF();
//	for (std::size_t i = 0; i < numClusters; ++i) {
//		const auto& cl1Reg = rNetwork.getCluster(i).getRegion();
//		Composition lo1 = cl1Reg.getOrigin();
//
//		std::cout << lo1[Spec::He] << ", " << lo1[Spec::D] << ", "
//				<< lo1[Spec::T] << ", " << lo1[Spec::V] << ", " << lo1[Spec::I]
//				<< std::endl;
//	}
//	std::cout << numClusters << std::endl;

	Composition comp;
	// Initialize the composition
	for (auto i : rNetwork.getSpeciesRange()) {
		comp[i] = 0;
	}
	comp[Spec::He] = 3;
	comp[Spec::V] = 2;
	comp[Spec::I] = 0;
	auto cluster = rNetwork.findCluster(comp, plsm::onHost);
	auto compRegion = cluster.getRegion();
	BOOST_REQUIRE(compRegion.isSimplex());
	BOOST_REQUIRE_EQUAL(compRegion.getOrigin(), comp);

//	cluster = rNetwork.getCluster(0, plsm::onHost);
//	compRegion = cluster.getRegion();
//	BOOST_REQUIRE_EQUAL(compRegion.getOrigin(), Composition { });

	using Reaction = typename NetworkType::ReactionType;
	// Reaction r0;
	// Reaction r(rNetwork, 0, Reaction::Type::dissociation, 0, 1, 2);
	// BOOST_REQUIRE(r.getType() == Reaction::Type::dissociation);
	auto concs = new double[3] { };
	auto fluxes = new double[3] { };
	using Vue = Reaction::ConcentrationsView;
	// r.contributeFlux(Vue(concs, 3), Vue(fluxes, 3), 0);
}

BOOST_AUTO_TEST_SUITE_END()
