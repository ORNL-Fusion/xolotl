#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>

#include <xolotl/core/network/AlloyTraits.h>
#include <xolotl/core/network/FeTraits.h>
#include <xolotl/core/network/NETraits.h>
#include <xolotl/core/network/PSITraits.h>
#include <xolotl/core/network/SpeciesEnumSequence.h>
#include <xolotl/core/network/ZrTraits.h>

using namespace std;

/**
 * This suite is responsible for testing the SpeciesEnumSequence.
 */
BOOST_AUTO_TEST_SUITE(SpeciesEnumSequence_testSuite)

/**
 * Method checking alloy.
 */
BOOST_AUTO_TEST_CASE(alloyCase)
{
	// Define the sequence
	using GroupingRange = xolotl::core::network::SpeciesForGrouping<
		xolotl::core::network::AlloySpecies, 6>;

	// It should always return 0 because they are all orthogonal
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::AlloySpecies::V),
		0);
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::AlloySpecies::Void),
		0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(
							xolotl::core::network::AlloySpecies::Faulted),
		0);
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::AlloySpecies::I),
		0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(
							xolotl::core::network::AlloySpecies::Perfect),
		0);
	BOOST_REQUIRE_EQUAL(GroupingRange::mapToMomentId(
							xolotl::core::network::AlloySpecies::Frank),
		0);

	return;
}

/**
 * Method checking iron.
 */
BOOST_AUTO_TEST_CASE(feCase)
{
	// Define the sequence
	using GroupingRange = xolotl::core::network::SpeciesForGrouping<
		xolotl::core::network::FeSpeciesList, 3>;

	// I goes on V
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::FeSpeciesList::He),
		0);
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::FeSpeciesList::V),
		1);
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::FeSpeciesList::I),
		1);

	return;
}

/**
 * Method checking PSI.
 */
BOOST_AUTO_TEST_CASE(psiCase)
{
	// Define the sequence
	using FullGroupingRange = xolotl::core::network::SpeciesForGrouping<
		xolotl::core::network::PSIFullSpeciesList, 5>;
	// I goes on V
	BOOST_REQUIRE_EQUAL(FullGroupingRange::mapToMomentId(
							xolotl::core::network::PSIFullSpeciesList::He),
		0);
	BOOST_REQUIRE_EQUAL(FullGroupingRange::mapToMomentId(
							xolotl::core::network::PSIFullSpeciesList::D),
		1);
	BOOST_REQUIRE_EQUAL(FullGroupingRange::mapToMomentId(
							xolotl::core::network::PSIFullSpeciesList::T),
		2);
	BOOST_REQUIRE_EQUAL(FullGroupingRange::mapToMomentId(
							xolotl::core::network::PSIFullSpeciesList::V),
		3);
	BOOST_REQUIRE_EQUAL(FullGroupingRange::mapToMomentId(
							xolotl::core::network::PSIFullSpeciesList::I),
		3);

	// Define the sequence
	using DeuteriumGroupingRange = xolotl::core::network::SpeciesForGrouping<
		xolotl::core::network::PSIDeuteriumSpeciesList, 4>;
	// I goes on V
	BOOST_REQUIRE_EQUAL(DeuteriumGroupingRange::mapToMomentId(
							xolotl::core::network::PSIDeuteriumSpeciesList::He),
		0);
	BOOST_REQUIRE_EQUAL(DeuteriumGroupingRange::mapToMomentId(
							xolotl::core::network::PSIDeuteriumSpeciesList::D),
		1);
	BOOST_REQUIRE_EQUAL(DeuteriumGroupingRange::mapToMomentId(
							xolotl::core::network::PSIDeuteriumSpeciesList::V),
		2);
	BOOST_REQUIRE_EQUAL(DeuteriumGroupingRange::mapToMomentId(
							xolotl::core::network::PSIDeuteriumSpeciesList::I),
		2);

	// Define the sequence
	using TritiumGroupingRange = xolotl::core::network::SpeciesForGrouping<
		xolotl::core::network::PSITritiumSpeciesList, 4>;
	// I goes on V
	BOOST_REQUIRE_EQUAL(TritiumGroupingRange::mapToMomentId(
							xolotl::core::network::PSITritiumSpeciesList::He),
		0);
	BOOST_REQUIRE_EQUAL(TritiumGroupingRange::mapToMomentId(
							xolotl::core::network::PSITritiumSpeciesList::T),
		1);
	BOOST_REQUIRE_EQUAL(TritiumGroupingRange::mapToMomentId(
							xolotl::core::network::PSITritiumSpeciesList::V),
		2);
	BOOST_REQUIRE_EQUAL(TritiumGroupingRange::mapToMomentId(
							xolotl::core::network::PSITritiumSpeciesList::I),
		2);

	// Define the sequence
	using HeliumGroupingRange = xolotl::core::network::SpeciesForGrouping<
		xolotl::core::network::PSIHeliumSpeciesList, 3>;
	// I goes on V
	BOOST_REQUIRE_EQUAL(HeliumGroupingRange::mapToMomentId(
							xolotl::core::network::PSIHeliumSpeciesList::He),
		0);
	BOOST_REQUIRE_EQUAL(HeliumGroupingRange::mapToMomentId(
							xolotl::core::network::PSIHeliumSpeciesList::V),
		1);
	BOOST_REQUIRE_EQUAL(HeliumGroupingRange::mapToMomentId(
							xolotl::core::network::PSIHeliumSpeciesList::I),
		1);

	return;
}

/**
 * Method checking alphaZr.
 */
BOOST_AUTO_TEST_CASE(zrCase)
{
	// Define the sequence
	using GroupingRange = xolotl::core::network::SpeciesForGrouping<
		xolotl::core::network::ZrSpecies, 3>;

	// It should always return 0 because they are all orthogonal
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::ZrSpecies::V), 0);
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::ZrSpecies::Basal),
		0);
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::ZrSpecies::I), 0);

	return;
}

/**
 * Method checking UO2.
 */
BOOST_AUTO_TEST_CASE(neCase)
{
	// Define the sequence
	using GroupingRange = xolotl::core::network::SpeciesForGrouping<
		xolotl::core::network::NESpecies, 1>;

	// It should always return 0 because they are all orthogonal
	BOOST_REQUIRE_EQUAL(
		GroupingRange::mapToMomentId(xolotl::core::network::NESpecies::Xe), 0);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
