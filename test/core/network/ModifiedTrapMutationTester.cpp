#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <algorithm>
#include <fstream>
#include <iostream>

#include <boost/test/unit_test.hpp>

#include <xolotl/core/network/PSIReactionNetwork.h>
#include <xolotl/core/network/impl/TrapMutationReaction.tpp>
#include <xolotl/test/CommandLine.h>
#include <xolotl/options/Options.h>

using namespace std;
using namespace xolotl;
using namespace core;

using Kokkos::ScopeGuard;
BOOST_GLOBAL_FIXTURE(ScopeGuard);

class TungstenTMTestHelper
{
public:
	using NetworkType =
		network::PSIReactionNetwork<network::PSIFullSpeciesList>;
	using AmountType = typename NetworkType::AmountType;
	using IndexType = NetworkType::IndexType;
	using TMReaction = typename NetworkType::Traits::TrapMutationReactionType;

	TungstenTMTestHelper(const std::string& materialName) :
		_grid(makeGrid(_nGrid)),
		_network(makeNetwork(materialName, _grid)),
		_dof(_network.getDOF() + 1),
		_nPartials(_network.getDiagonalFill(_dfill)),
		_indices(_dof)
	{
		for (IndexType i = 0, id = 0; i < _dof; ++i) {
			const auto& row = _dfill[i];
			for (IndexType j = 0; j < row.size(); ++j) {
				_indices[i].push_back(id);
				++id;
			}
		}
	}

	std::pair<double, double>
	getDepthAndSpacing(int i) const
	{
		auto surfacePos = _grid[_surfacePosition + 1];
		if (i == _grid.size() - 1) {
			auto depth = _grid[i] - surfacePos;
			auto spacing = _grid[i] - _grid[i - 1];
			return std::make_pair(depth, spacing);
		}
		auto curXPos = (_grid[i] + _grid[i + 1]) / 2.0;
		auto prevXPos = (_grid[i - 1] + _grid[i]) / 2.0;
		auto depth = curXPos - surfacePos;
		auto spacing = curXPos - prevXPos;
		return std::make_pair(depth, spacing);
	}

	bool
	dfillFind(int row, int col) const
	{
		auto& dfillRow = _dfill.at(row);
		auto it = std::find(begin(dfillRow), end(dfillRow), col);
		return it != end(dfillRow);
	}

	IndexType
	getPartialsIndex(int row, int col) const
	{
		auto& dfillRow = _dfill.at(row);
		auto it = std::find(begin(dfillRow), end(dfillRow), col);
		return _indices[row][std::distance(begin(dfillRow), it)];
	}

	void
	setTemperatures(double temp)
	{
		std::vector<double> temperatures(_nGrid, temp);
		_network.setTemperatures(temperatures);
		_network.syncClusterDataOnHost();
	}

	Kokkos::View<double*>
	getConcentrations(int gridIndex) const
	{
		std::vector<double> concentrations(_dof);
		for (std::size_t i = 0; i < _dof; ++i) {
			double fac = static_cast<double>(gridIndex) * _dof + i;
			concentrations[i] = fac * fac;
		}
		using HostUnmanaged =
			Kokkos::View<double*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
		auto hConcs = HostUnmanaged(concentrations.data(), _dof);
		auto dConcs = Kokkos::View<double*>("Concentrations", _dof);
		deep_copy(dConcs, hConcs);
		return dConcs;
	}

	Kokkos::View<double*, Kokkos::HostSpace>
	computeFluxes(int gridIndex)
	{
		auto dConcs = getConcentrations(gridIndex);
		auto dFluxes = Kokkos::View<double*>("Fluxes", _dof);
		auto [curDepth, curSpacing] = getDepthAndSpacing(gridIndex);
		_network.computeFluxes<TMReaction>(
			dConcs, dFluxes, gridIndex, curDepth, curSpacing);
		auto hFluxes = create_mirror_view(dFluxes);
		deep_copy(hFluxes, dFluxes);

		return hFluxes;
	}

	Kokkos::View<double*, Kokkos::HostSpace>
	computePartials(int gridIndex)
	{
		auto dConcs = getConcentrations(gridIndex);
		auto dPartials = Kokkos::View<double*>("Partials", _nPartials);
		auto [curDepth, curSpacing] = getDepthAndSpacing(gridIndex);
		_network.computePartials<TMReaction>(
			dConcs, dPartials, gridIndex, curDepth, curSpacing);
		auto hPartials = create_mirror_view(dPartials);
		deep_copy(hPartials, dPartials);

		return hPartials;
	}

private:
	static std::vector<double>
	makeGrid(int nGrid)
	{
		// Suppose we have a grid with 13 grid points and distance of
		// 0.1 nm between grid points
		std::vector<double> grid;
		for (int l = 0; l < nGrid; l++) {
			grid.push_back((double)l * 0.1);
		}

		return grid;
	}

	static NetworkType
	makeNetwork(
		const std::string& materialName, const std::vector<double>& grid)
	{
		// Create the option to create a network
		xolotl::options::Options opts;
		// Create a good parameter file
		std::string parameterFile = "param.txt";
		std::ofstream paramFile(parameterFile);
		paramFile << "netParam=8 0 0 10 6" << '\n'
				  << "process=reaction modifiedTM" << '\n'
				  << "material=" << materialName << '\n';
		paramFile.close();

		// Create a fake command line to read the options
        test::CommandLine<2> cl{{"fakeXolotlAppNameForTests", parameterFile}};
        opts.readParams(cl.argc, cl.argv);

		// Remove the created file
		std::remove(parameterFile.c_str());

		AmountType maxV = opts.getMaxV();
		AmountType maxI = opts.getMaxI();
		AmountType maxHe = opts.getMaxImpurity();
		AmountType maxD = opts.getMaxD();
		AmountType maxT = opts.getMaxT();
		NetworkType network({maxHe, maxD, maxT, maxV, maxI}, grid.size(), opts);
		network.syncClusterDataOnHost();
		network.getSubpaving().syncZones(plsm::onHost);

		return network;
	}

private:
	int _nGrid{13};
	int _surfacePosition{0};
	std::vector<double> _grid;
	NetworkType _network;
	int _dof{};
	NetworkType::SparseFillMap _dfill;
	IndexType _nPartials;
	std::vector<std::vector<IndexType>> _indices;
};

/**
 * This suite is responsible for testing the W100TrapMutationHandler.
 */
BOOST_AUTO_TEST_SUITE(mtm)

/**
 * Method checking the initialization and the compute modified trap-mutation
 * methods.
 */
BOOST_AUTO_TEST_CASE(W100)
{
	TungstenTMTestHelper hlp("W100");

	// Check some values in dfill
	BOOST_REQUIRE(hlp.dfillFind(27, 27));
	BOOST_REQUIRE(hlp.dfillFind(0, 27));
	BOOST_REQUIRE(hlp.dfillFind(28, 27));
	BOOST_REQUIRE(hlp.dfillFind(38, 38));
	BOOST_REQUIRE(hlp.dfillFind(39, 38));
	BOOST_REQUIRE(hlp.dfillFind(0, 38));
	BOOST_REQUIRE(hlp.dfillFind(60, 60));
	BOOST_REQUIRE(hlp.dfillFind(61, 60));
	BOOST_REQUIRE(hlp.dfillFind(0, 60));
	BOOST_REQUIRE(hlp.dfillFind(71, 71));
	BOOST_REQUIRE(hlp.dfillFind(73, 71));
	BOOST_REQUIRE(hlp.dfillFind(1, 71));

	// Set the temperature to compute the rates
	hlp.setTemperatures(1200.0);

	// Compute the modified trap mutation at the sixth grid point
	auto fluxes = hlp.computeFluxes(6);
	BOOST_REQUIRE_CLOSE(fluxes[0], 2.25691e+24, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(fluxes[27], -2.25691e+24, 0.01); // He2
	BOOST_REQUIRE_CLOSE(fluxes[28], 2.25691e+24, 0.01); // Create He2V

	// Compute the modified trap mutation at the ninth grid point
	fluxes = hlp.computeFluxes(9);
	BOOST_REQUIRE_CLOSE(fluxes[1], 1.3724e+21, 0.01); // Create I2
	BOOST_REQUIRE_CLOSE(fluxes[38], 0.0, 0.01); // He3
	BOOST_REQUIRE_CLOSE(fluxes[39], 0.0, 0.01); // Doesn't create He3V
	BOOST_REQUIRE_CLOSE(fluxes[82], -6.9358e+20, 0.01); // He7
	BOOST_REQUIRE_CLOSE(fluxes[84], 6.9358e+20, 0.01); // Create He7V2

	// Compute the partials for the modified trap-mutation at the 9th grid point
	auto partials = hlp.computePartials(9);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(82, 82)], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(84, 82)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(1, 82)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(71, 71)], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(73, 71)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(1, 71)], 6.575931697e+14, 0.01);

	// Change the temperature of the network
	hlp.setTemperatures(500.0);

	// Compute the partials for the modified trap-mutation at the 9th grid point
	partials = hlp.computePartials(9);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(60, 60)], -5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(61, 60)], 5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 60)], 5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(71, 71)], -5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(73, 71)], 5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(1, 71)], 5.53624e+14, 0.01);
}

BOOST_AUTO_TEST_CASE(W110)
{
	TungstenTMTestHelper hlp("W110");

	// Check some values in dfill
	BOOST_REQUIRE(hlp.dfillFind(27, 27));
	BOOST_REQUIRE(hlp.dfillFind(0, 27));
	BOOST_REQUIRE(hlp.dfillFind(28, 27));
	BOOST_REQUIRE(hlp.dfillFind(38, 38));
	BOOST_REQUIRE(hlp.dfillFind(39, 38));
	BOOST_REQUIRE(hlp.dfillFind(0, 38));
	BOOST_REQUIRE(hlp.dfillFind(60, 60));
	BOOST_REQUIRE(hlp.dfillFind(61, 60));
	BOOST_REQUIRE(hlp.dfillFind(0, 60));
	BOOST_REQUIRE(hlp.dfillFind(71, 71));
	BOOST_REQUIRE(hlp.dfillFind(73, 71));
	BOOST_REQUIRE(hlp.dfillFind(1, 71));
	BOOST_REQUIRE(hlp.dfillFind(49, 49));
	BOOST_REQUIRE(hlp.dfillFind(50, 49));
	BOOST_REQUIRE(hlp.dfillFind(0, 49));

	// Set the temperature to compute the rates
	hlp.setTemperatures(1200.0);

	// Compute the modified trap mutation at the eigth grid point
	auto fluxes = hlp.computeFluxes(8);
	BOOST_REQUIRE_CLOSE(fluxes[0], 3.370929e+24, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(fluxes[27], -3.370929e+24, 0.01); // He2
	BOOST_REQUIRE_CLOSE(fluxes[28], 3.370929e+24, 0.01); // Create He2V

	// Compute the modified trap mutation at the tenth grid point
	fluxes = hlp.computeFluxes(10);
	BOOST_REQUIRE_CLOSE(fluxes[0], 2.38288e+21, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(fluxes[27], 0.0, 0.01); // He2
	BOOST_REQUIRE_CLOSE(fluxes[28], 0.0, 0.01); // Doesn't create He2V
	BOOST_REQUIRE_CLOSE(fluxes[60], -8.1022e+20, 0.01); // He5
	BOOST_REQUIRE_CLOSE(fluxes[61], 8.1022e+20, 0.01); // Create He5V

	// Compute the partials for the modified trap-mutation at the 9th grid point
	auto partials = hlp.computePartials(10);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(38, 38)], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(39, 38)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 38)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(49, 49)], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(50, 49)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 49)], 6.575931697e+14, 0.01);

	// Change the temperature of the network
	hlp.setTemperatures(500.0);

	// Compute the partials for the modified trap-mutation at the 9th grid point
	partials = hlp.computePartials(10);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(38, 38)], -5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(39, 38)], 5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 38)], 5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(49, 49)], -5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(50, 49)], 5.53624e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 49)], 5.53624e+14, 0.01);
}

BOOST_AUTO_TEST_CASE(W111)
{
	TungstenTMTestHelper hlp("W111");

	// Check some values in dfill
	BOOST_REQUIRE(hlp.dfillFind(27, 27));
	BOOST_REQUIRE(hlp.dfillFind(0, 27));
	BOOST_REQUIRE(hlp.dfillFind(28, 27));
	BOOST_REQUIRE(hlp.dfillFind(16, 16));
	BOOST_REQUIRE(hlp.dfillFind(17, 16));
	BOOST_REQUIRE(hlp.dfillFind(0, 16));
	BOOST_REQUIRE(hlp.dfillFind(38, 38));
	BOOST_REQUIRE(hlp.dfillFind(39, 38));
	BOOST_REQUIRE(hlp.dfillFind(0, 38));
	BOOST_REQUIRE(hlp.dfillFind(49, 49));
	BOOST_REQUIRE(hlp.dfillFind(50, 49));
	BOOST_REQUIRE(hlp.dfillFind(0, 49));

	// Set the temperature to compute the rates
	hlp.setTemperatures(1200.0);

	// Compute the modified trap mutation at the seventh grid point
	auto fluxes = hlp.computeFluxes(7);
	BOOST_REQUIRE_CLOSE(fluxes[0], 3.315351e+24, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(fluxes[16], -3.315351e+24, 0.01); // He2
	BOOST_REQUIRE_CLOSE(fluxes[17], 3.315351e+24, 0.01); // Create He2V

	// Compute the modified trap mutation at the twelfth grid point
	fluxes = hlp.computeFluxes(12);
	BOOST_REQUIRE_CLOSE(fluxes[0], 4.545446e+21, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(fluxes[27], 0.0, 0.01); // He2
	BOOST_REQUIRE_CLOSE(fluxes[28], 0.0, 0.01); // Doesn't create He2V
	BOOST_REQUIRE_CLOSE(fluxes[49], -1.12677e+21, 0.01); // He4
	BOOST_REQUIRE_CLOSE(fluxes[50], 1.12677e+21, 0.01); // Create He4V

	// Compute the partials for the modified trap-mutation at the twelfth grid
	// point
	auto partials = hlp.computePartials(12);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(38, 38)], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(39, 38)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 38)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(49, 49)], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(50, 49)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 49)], 6.575931697e+14, 0.01);

	// Change the temperature of the network
	hlp.setTemperatures(500.0);

	// Compute the partials for the modified trap-mutation at the twelfth grid
	// point
	partials = hlp.computePartials(12);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(38, 38)], -5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(39, 38)], 5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 38)], 5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(49, 49)], -5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(50, 49)], 5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 49)], 5.536237e+14, 0.01);
}

BOOST_AUTO_TEST_CASE(W211)
{
	TungstenTMTestHelper hlp("W211");

	// Check some values in dfill
	BOOST_REQUIRE(hlp.dfillFind(27, 27));
	BOOST_REQUIRE(hlp.dfillFind(0, 27));
	BOOST_REQUIRE(hlp.dfillFind(28, 27));
	BOOST_REQUIRE(hlp.dfillFind(16, 16));
	BOOST_REQUIRE(hlp.dfillFind(17, 16));
	BOOST_REQUIRE(hlp.dfillFind(0, 16));
	BOOST_REQUIRE(hlp.dfillFind(38, 38));
	BOOST_REQUIRE(hlp.dfillFind(39, 38));
	BOOST_REQUIRE(hlp.dfillFind(0, 38));
	BOOST_REQUIRE(hlp.dfillFind(49, 49));
	BOOST_REQUIRE(hlp.dfillFind(51, 49));
	BOOST_REQUIRE(hlp.dfillFind(1, 49));

	// Set the temperature to compute the rates
	hlp.setTemperatures(1200.0);

	// Compute the modified trap mutation at the sixth grid point
	auto fluxes = hlp.computeFluxes(6);
	BOOST_REQUIRE_CLOSE(fluxes[0], 6.88681e+23, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(fluxes[16], -6.88681e+23, 0.01); // He
	BOOST_REQUIRE_CLOSE(fluxes[17], 6.88681e+23, 0.01); // Create HeV

	// Compute the modified trap mutation at the eleventh grid point
	fluxes = hlp.computeFluxes(11);
	BOOST_REQUIRE_CLOSE(fluxes[0], 9.359188e+20, 0.01); // Create I
	BOOST_REQUIRE_CLOSE(fluxes[27], 0.0, 0.01); // He2
	BOOST_REQUIRE_CLOSE(fluxes[28], 0.0, 0.01); // Doesn't create He2V
	BOOST_REQUIRE_CLOSE(fluxes[49], -9.53258e+20, 0.01); // He4
	BOOST_REQUIRE_CLOSE(fluxes[51], 9.53258e+20, 0.01); // Create He4V2

	// Compute the partials for the modified trap-mutation at the eleventh grid
	// point
	auto partials = hlp.computePartials(11);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(38, 38)], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(39, 38)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 38)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(49, 49)], -6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(51, 49)], 6.575931697e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(1, 49)], 6.575931697e+14, 0.01);

	// Change the temperature of the network
	hlp.setTemperatures(500.0);

	// Compute the partials for the modified trap-mutation at the eleventh grid
	// point
	partials = hlp.computePartials(11);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(38, 38)], -5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(39, 38)], 5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(0, 38)], 5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(49, 49)], -5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(51, 49)], 5.536237e+14, 0.01);
	BOOST_REQUIRE_CLOSE(
		partials[hlp.getPartialsIndex(1, 49)], 5.536237e+14, 0.01);
}

BOOST_AUTO_TEST_SUITE_END()
