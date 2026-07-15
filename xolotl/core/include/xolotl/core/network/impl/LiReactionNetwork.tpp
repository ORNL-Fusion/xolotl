#pragma once

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/detail/impl/TrapReactionGenerator.tpp>
#include <xolotl/core/network/impl/LiClusterGenerator.tpp>
#include <xolotl/core/network/impl/LiReaction.tpp>
#include <xolotl/core/network/impl/ReactionNetwork.tpp>
#include <xolotl/util/Tokenizer.h>

namespace xolotl
{
namespace core
{
namespace network
{
double
LiReactionNetwork::computeTritiumRelease(
	double thetaH, double thetaT, double TSphereConc, double temperature)
{
	double RT = getRT(temperature);

	// Get the H concentration in the purge
	auto Hconc = computeHydrogenPurge(temperature);

	// Get the T concentration in the purge
	auto Tconc = computeTritiumPurge(thetaH, thetaT, Hconc, temperature);

	// Get the molecule concentrations in the purge
	auto T2conc = Tconc * Tconc / (2.0 * (Hconc + Tconc));
	auto HTconc = Hconc * Tconc / (Hconc + Tconc);

	// Compute the desorption rate
	auto qads = getQads(thetaH + thetaT);
	double desorptionE = qads + adsorptionE;
	double desorptionRate =
		(surfaceSites * adjSites * RT * exp(-2.0 * desorptionE / RT)) /
		(2.0 * nAv * planckConstant);

	// Compute the adsorption rate
	double adsorptionRate = (1.0 - thetaH - thetaT) * (1.0 - thetaH - thetaT) *
		condensationCoef * adjSites * exp(-2.0 * adsorptionE / RT) *
		(gamma_HT_T * HTconc * sqrt(RT / htMolWeight) +
			gamma_T2_T * T2conc * sqrt(RT / t2MolWeight)) /
		sqrt(8.0e-3 * ::xolotl::core::pi);

	return sphereSurface *
		(desorptionRate * (gamma_HT_T * thetaH + gamma_T2_T * thetaT) * thetaT -
			adsorptionRate);
}

double
LiReactionNetwork::getTritiumFlux(
	double thetaH, double thetaT, double TSphereConc, double temperature)
{
	double RT = getRT(temperature);

	// Beta rate
	double tFlux = TSphereConc * (1.0 - thetaH - thetaT) * 1.0e13 *
		exp(-betaE / RT) / sqrt(surfaceSites);

	// Dissolution rate
	auto qads = getQads(thetaH + thetaT);
	double dissolutionE = solutionE + betaE + qads;
	tFlux -= surfaceSites * adjSites * RT * exp(-dissolutionE / RT) * thetaT /
		(2.0 * nAv * planckConstant);

	return -tFlux;
}

void
LiReactionNetwork::initializeExtraDOFs(const options::IOptions& options)
{
	this->_clusterData.h_view().setHCoverageId(this->_numDOFs);
	this->_clusterData.h_view().setTCoverageId(this->_numDOFs + 1);
	this->_numDOFs += 2;
}

void
LiReactionNetwork::computeFluxesPreProcess(ConcentrationsView concentrations,
	FluxesView fluxes, IndexType gridIndex, double surfaceDepth, double spacing)
{
	if (util::equal(surfaceDepth, totalDepth)) {
		// Get the temperature
		auto temperature = this->_clusterData.h_view().temperature(gridIndex);
		double RT = getRT(temperature);

		// Get the H concentration in the purge
		auto Hconc = computeHydrogenPurge(temperature);

		// Get the surface coverages
		auto thetaH = concentrations(this->_clusterData.h_view().hCoverageId());
		auto thetaT = concentrations(this->_clusterData.h_view().tCoverageId());

		// Get the T concentration in the purge
		auto Tconc = computeTritiumPurge(thetaH, thetaT, Hconc, temperature);

		// Get the molecule concentrations in the purge
		auto T2conc = Tconc * Tconc / (2.0 * (Hconc + Tconc));
		auto H2conc = Hconc * Hconc / (2.0 * (Hconc + Tconc));
		auto HTconc = Hconc * Tconc / (Hconc + Tconc);

		// Rate equation for thetaH
		// Adsorption term
		double hFlux = (1.0 - thetaH - thetaT) * (1.0 - thetaH - thetaT) *
			condensationCoef * adjSites * exp(-2.0 * adsorptionE / RT) *
			(gamma_HT_H * HTconc * sqrt(RT / htMolWeight) +
				gamma_H2_H * H2conc * sqrt(RT / h2MolWeight)) /
			sqrt(8.0e-3 * ::xolotl::core::pi);

		// Desorption term
		auto qads = getQads(thetaH + thetaT);
		double desorptionE = qads + adsorptionE;
		double desorptionRate =
			(surfaceSites * adjSites * RT * exp(-2.0 * desorptionE / RT)) /
			(2.0 * nAv * planckConstant);
		hFlux -= desorptionRate *
			(gamma_HT_H * thetaH * thetaT + gamma_H2_H * thetaH * thetaH);

		// Rate equation for thetaT
		// Adsorption term
		double tFlux = (1.0 - thetaH - thetaT) * (1.0 - thetaH - thetaT) *
			condensationCoef * adjSites * exp(-2.0 * adsorptionE / RT) *
			(gamma_HT_T * HTconc * sqrt(RT / htMolWeight) +
				gamma_T2_T * T2conc * sqrt(RT / t2MolWeight)) /
			sqrt(8.0e-3 * ::xolotl::core::pi);

		// Desorption term
		tFlux -= desorptionRate *
			(gamma_HT_T * thetaH * thetaT + gamma_T2_T * thetaT * thetaT);

		// Concentration of T at the edge of the sphere
		auto TSphereConc = concentrations(0);
		// Beta term
		tFlux += TSphereConc * (1.0 - thetaH - thetaT) * 1.0e13 *
			exp(-betaE / RT) / sqrt(surfaceSites);

		// Dissolution term
		double dissolutionE = solutionE + betaE + qads;
		tFlux -= surfaceSites * adjSites * RT * exp(-dissolutionE / RT) *
			thetaT / (2.0 * nAv * planckConstant);

		// Update the RHS function
		Kokkos::atomic_add(&fluxes[this->_clusterData.h_view().hCoverageId()],
			hFlux / surfaceSites);
		Kokkos::atomic_add(&fluxes[this->_clusterData.h_view().tCoverageId()],
			tFlux / surfaceSites);
	}
}

void
LiReactionNetwork::computePartialsPreProcess(ConcentrationsView concentrations,
	Kokkos::View<double*> values, IndexType gridIndex, double surfaceDepth,
	double spacing)
{
	if (util::equal(surfaceDepth, totalDepth)) {
		// Get the temperature
		auto temperature = this->_clusterData.h_view().temperature(gridIndex);
		double RT = getRT(temperature);

		// Get the H concentration in the purge
		auto Hconc = computeHydrogenPurge(temperature);

		// Get the surface coverages
		auto thetaH = concentrations(this->_clusterData.h_view().hCoverageId());
		auto thetaT = concentrations(this->_clusterData.h_view().tCoverageId());

		// Get the T concentration in the purge
		auto Tconc = computeTritiumPurge(thetaH, thetaT, Hconc, temperature);
		auto TconcThetaH =
			computeTritiumPurgeThetaH(thetaH, thetaT, Hconc, temperature);
		auto TconcThetaT =
			computeTritiumPurgeThetaT(thetaH, thetaT, Hconc, temperature);

		// Get the molecule concentrations in the purge
		auto T2conc = Tconc * Tconc / (2.0 * (Hconc + Tconc));
		auto H2conc = Hconc * Hconc / (2.0 * (Hconc + Tconc));
		auto HTconc = Hconc * Tconc / (Hconc + Tconc);
		auto T2concThetaH = (2.0 * TconcThetaH * Tconc * (Hconc + Tconc) -
								TconcThetaH * Tconc * Tconc) /
			(2.0 * (Hconc + Tconc) * (Hconc + Tconc));
		auto H2concThetaH = -(Hconc * Hconc * TconcThetaH) /
			(2.0 * (Hconc + Tconc) * (Hconc + Tconc));
		auto HTconcThetaH = (Hconc * TconcThetaH * (Hconc + Tconc) -
								Hconc * Tconc * TconcThetaH) /
			((Hconc + Tconc) * (Hconc + Tconc));
		auto T2concThetaT = (2.0 * TconcThetaT * Tconc * (Hconc + Tconc) -
								TconcThetaT * Tconc * Tconc) /
			(2.0 * (Hconc + Tconc) * (Hconc + Tconc));
		auto H2concThetaT = -(Hconc * Hconc * TconcThetaT) /
			(2.0 * (Hconc + Tconc) * (Hconc + Tconc));
		auto HTconcThetaT = (Hconc * TconcThetaT * (Hconc + Tconc) -
								Hconc * Tconc * TconcThetaT) /
			((Hconc + Tconc) * (Hconc + Tconc));

		// Rate equation for thetaH
		// Adsorption term
		double alpha = condensationCoef * adjSites *
			exp(-2.0 * adsorptionE / RT) / sqrt(8.0e-3 * ::xolotl::core::pi);
		double beta = gamma_HT_H * sqrt(RT / htMolWeight);
		double gamma = gamma_H2_H * sqrt(RT / h2MolWeight);
		double hPartialThetaH = -2.0 * alpha * (1.0 - thetaH - thetaT) *
				(HTconc * beta + H2conc * gamma) +
			(1.0 - thetaH - thetaT) * (1.0 - thetaH - thetaT) * alpha *
				(beta * HTconcThetaH + gamma * H2concThetaH);
		double hPartialThetaT = -2.0 * alpha * (1.0 - thetaH - thetaT) *
				(HTconc * beta + H2conc * gamma) +
			(1.0 - thetaH - thetaT) * (1.0 - thetaH - thetaT) * alpha *
				(beta * HTconcThetaT + gamma * H2concThetaT);

		// Desorption term
		auto qads = getQads(thetaH + thetaT);
		double desorptionE = qads + adsorptionE;
		double desorptionRate =
			(surfaceSites * adjSites * RT * exp(-2.0 * desorptionE / RT)) /
			(2.0 * nAv * planckConstant);
		double desorptionRateTheta = -surfaceSites * adjSites *
			getQadsTheta(thetaH + thetaT) * exp(-2.0 * desorptionE / RT) /
			(nAv * planckConstant);
		hPartialThetaH -= desorptionRateTheta *
				(gamma_HT_H * thetaH * thetaT + gamma_H2_H * thetaH * thetaH) +
			desorptionRate * (gamma_HT_H * thetaT + 2.0 * gamma_H2_H * thetaH);
		hPartialThetaT -= desorptionRateTheta *
				(gamma_HT_H * thetaH * thetaT + gamma_H2_H * thetaH * thetaH) +
			desorptionRate * gamma_HT_H * thetaH;

		// Update the partial derivative of thetaH with thetaH and thetaT
		Kokkos::atomic_add(
			&values[this->_connEntries(0)], hPartialThetaH / surfaceSites);
		Kokkos::atomic_add(
			&values[this->_connEntries(1)], hPartialThetaT / surfaceSites);

		// Rate equation for thetaT
		// Adsorption term
		beta = gamma_HT_T * sqrt(RT / htMolWeight);
		double delta = gamma_T2_T * sqrt(RT / t2MolWeight);
		double tPartialThetaH = -2.0 * alpha * (1.0 - thetaH - thetaT) *
				(HTconc * beta + T2conc * delta) +
			(1.0 - thetaH - thetaT) * (1.0 - thetaH - thetaT) * alpha *
				(beta * HTconcThetaH + delta * T2concThetaH);
		double tPartialThetaT = -2.0 * alpha * (1.0 - thetaH - thetaT) *
				(HTconc * beta + T2conc * delta) +
			(1.0 - thetaH - thetaT) * (1.0 - thetaH - thetaT) * alpha *
				(beta * HTconcThetaT + delta * T2concThetaT);

		// Desorption term
		tPartialThetaH -= desorptionRateTheta *
				(gamma_HT_T * thetaH * thetaT + gamma_T2_T * thetaT * thetaT) +
			desorptionRate * gamma_HT_T * thetaT;
		tPartialThetaT -= desorptionRateTheta *
				(gamma_HT_T * thetaH * thetaT + gamma_T2_T * thetaT * thetaT) +
			desorptionRate * (gamma_HT_T * thetaH + 2.0 * gamma_T2_T * thetaT);

		// Concentration of T at the edge of the sphere
		auto TSphereConc = concentrations(0);
		// Beta term
		tPartialThetaH -=
			TSphereConc * 1.0e13 * exp(-betaE / RT) / sqrt(surfaceSites);
		tPartialThetaT -=
			TSphereConc * 1.0e13 * exp(-betaE / RT) / sqrt(surfaceSites);

		// Dissolution term
		double dissolutionE = solutionE + betaE + qads;
		tPartialThetaH += surfaceSites * adjSites * exp(-dissolutionE / RT) *
			thetaT * getQadsTheta(thetaH + thetaT) /
			(2.0 * nAv * planckConstant);
		tPartialThetaT -= surfaceSites * adjSites * RT *
			exp(-dissolutionE / RT) *
			(1.0 - thetaT * getQadsTheta(thetaH + thetaT) / RT) /
			(2.0 * nAv * planckConstant);

		// Update the partial derivative of thetaT with thetaH and thetaT
		Kokkos::atomic_add(
			&values[this->_connEntries(2)], tPartialThetaH / surfaceSites);
		Kokkos::atomic_add(
			&values[this->_connEntries(3)], tPartialThetaT / surfaceSites);

		// Last term for sphere T concentration
		double tPartialConc = (1.0 - thetaH - thetaT) * 1.0e13 *
			exp(-betaE / RT) / sqrt(surfaceSites);
		Kokkos::atomic_add(
			&values[this->_connEntries(5)], tPartialConc / surfaceSites);
	}
}

KOKKOS_INLINE_FUNCTION
void
LiReactionNetwork::setConnectivity(LiReactionNetwork::Connectivity conn)
{
	// Count
	IndexType nPartials =
		4; // For the surface coverages, both interact with both

	auto& subpaving = this->getSubpaving();

	// Connects to single hydrogen and vice versa
	Composition comp = Composition::zero();
	comp[Species::H] = 1;
	auto clusterId = subpaving.findTileId(comp);
	if (clusterId != subpaving.invalidIndex())
		nPartials += 4;

	// Create the view
	this->_connEntries = Kokkos::View<IndexType*>(
		"Connectivity entries for extra DOFs", nPartials);

	// Save the entries
	auto hCovId = this->_clusterData.h_view().hCoverageId();
	auto tCovId = this->_clusterData.h_view().tCoverageId();

	// Surface coverages connect to themselves
	this->_connEntries(0) = conn(hCovId, hCovId);
	this->_connEntries(1) = conn(hCovId, tCovId);
	this->_connEntries(2) = conn(tCovId, hCovId);
	this->_connEntries(3) = conn(tCovId, tCovId);

	// Hydrogen connects to them
	this->_connEntries(4) = conn(hCovId, clusterId);
	this->_connEntries(5) = conn(tCovId, clusterId);
	this->_connEntries(6) = conn(clusterId, hCovId);
	this->_connEntries(7) = conn(clusterId, tCovId);

	return;
}

void
LiReactionNetwork::setReactionParams(
	std::vector<double> grid, std::string trapParamFile)
{
	// Read the parameter file
	std::ifstream paramFile(trapParamFile);

	if (!paramFile.good()) {
		// Print a message
		XOLOTL_LOG
			<< "No parameter file for trap reactions, they will not be used.";
	}
	else {
		// Get the reactions
		auto trapReactions =
			this->_reactions.template getView<LiTrapReaction>();

		// Get the line
		std::string line;
		getline(paramFile, line);

		// Read the first line
		auto tokens = util::Tokenizer<double>{line}();
		// And start looping on the lines
		int index = 0;
		while (tokens.size() > 0) {
			if (tokens.size() != 3) {
				throw std::runtime_error(
					"\nNot the correct number of trap parameters for "
					"the TrapReaction: 3 expected.");
			}

			// Keep the values
			auto params = tokens;
			std::vector<double> densities;

			// Get the parameters for the fit
			getline(paramFile, line);
			tokens = util::Tokenizer<double>{line}();
			if (tokens.size() == 1) {
				// Same density everywhere
				densities = std::vector<double>(this->_gridSize, tokens[0]);
			}
			else if (tokens.size() != 17) {
				throw std::runtime_error(
					"\nNot the correct number of fit parameters for the "
					"trap reaction density: 17 expected.");
			}
			else {
				double totalDepth = tokens[16] + 0.1;

				// Loop on the grid to provide the density at each grid point
				densities.push_back(0.0);
				for (auto j = 1; j < grid.size(); j++) {
					double den = 0.0;
					double x = (grid[j] + grid[j + 1]) / 2.0 - grid[1];
					if (x > totalDepth)
						den = 0.0;
					else {
						for (auto i = 0; i < 16; i++)
							den += tokens[i] * pow(x, (double)i);
					}

					densities.push_back(den);
				}
			}

			// Set the reaction parameters
			trapReactions(index).setParameters(
				densities, params[0], params[1], params[2]);
			index++;

			// Read the next line
			getline(paramFile, line);
			tokens = util::Tokenizer<double>{line}();
		}
	}
}

namespace detail
{
template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
LiReactionGenerator::operator()(IndexType i, IndexType j, TTag tag) const
{
	// Check the diffusion factors
	auto diffusionFactor = this->_clusterData.diffusionFactor;
	if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
		return;
	}

	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;
	using AmountType = typename NetworkType::AmountType;

	constexpr auto species = NetworkType::getSpeciesRange();
	constexpr auto speciesNoI = NetworkType::getSpeciesRangeNoI();

	// Add the trapping reactions between single mobile gas and traps
	addTraps(i, j, tag);

	auto numClusters = this->getNumberOfClusters();

	// Get the composition of each cluster
	const auto& cl1Reg = this->getCluster(i).getRegion();
	const auto& cl2Reg = this->getCluster(j).getRegion();
	Composition lo1 = cl1Reg.getOrigin();
	Composition hi1 = cl1Reg.getUpperLimitPoint();
	Composition lo2 = cl2Reg.getOrigin();
	Composition hi2 = cl2Reg.getUpperLimitPoint();

	auto& subpaving = this->getSubpaving();

	// General case
	constexpr auto numSpeciesNoI = NetworkType::getNumberOfSpeciesNoI();
	using BoundsArray =
		Kokkos::Array<Kokkos::pair<AmountType, AmountType>, numSpeciesNoI>;
	plsm::EnumIndexed<BoundsArray, Species> bounds;
	// Loop on the species
	for (auto l : species) {
		auto low = lo1[l] + lo2[l];
		auto high = hi1[l] + hi2[l] - 2;

		bounds[l] = {low, high};
	}

	// Look for potential product
	IndexType nProd = 0;
	for (IndexType k = 0; k < numClusters; ++k) {
		// Get the composition
		const auto& prodReg = this->getCluster(k).getRegion();
		bool isGood = true;
		// Loop on the species
		for (auto l : speciesNoI) {
			if (prodReg[l()].begin() > bounds[l()].second) {
				isGood = false;
				break;
			}
			if (prodReg[l()].end() - 1 < bounds[l()].first) {
				isGood = false;
				break;
			}
		}

		if (isGood) {
			// Increase nProd
			nProd++;
			this->addProductionReaction(tag, {i, j, k});
			// Loop on the species
			bool isOnAxis1 = false, isOnAxis2 = false;
			for (auto l : species) {
				if (lo1.isOnAxis(l()) && lo1[l()] == 1)
					isOnAxis1 = true;
				if (lo2.isOnAxis(l()) && lo2[l()] == 1)
					isOnAxis2 = true;
			}
			if (isOnAxis1 || isOnAxis2) {
				this->addDissociationReaction(tag, {k, i, j});
			}
		}
	}
}

void
LiReactionGenerator::addConnectivity(Connectivity& conn)
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	Kokkos::parallel_for(
		"PSIReactionGenerator::addConnectivity", 1, KOKKOS_LAMBDA(IndexType i) {
			auto& subpaving = this->getSubpaving();

			// Add connectivities for the surface coverage
			// Get the Ids
			auto hCovId = this->_clusterData.hCoverageId();
			auto tCovId = this->_clusterData.tCoverageId();
			// Connects to itself
			conn.add(hCovId, hCovId);
			conn.add(tCovId, tCovId);

			// Connects to each other
			conn.add(hCovId, tCovId);
			conn.add(tCovId, hCovId);

			// Connects to single hydrogen
			Composition comp = Composition::zero();
			comp[Species::H] = 1;
			auto clusterId = subpaving.findTileId(comp);
			if (clusterId != subpaving.invalidIndex()) {
				conn.add(hCovId, clusterId);
				conn.add(tCovId, clusterId);
				conn.add(clusterId, hCovId);
				conn.add(clusterId, tCovId);
			}
		});
}

template <typename TTag>
KOKKOS_INLINE_FUNCTION
void
LiReactionGenerator::addTraps(IndexType i, IndexType j, TTag tag) const
{
	using Species = typename NetworkType::Species;
	using Composition = typename NetworkType::Composition;

	// Get the compositions
	const auto& clReg1 = this->getCluster(i).getRegion();
	Composition lo1 = clReg1.getOrigin();
	const auto& clReg2 = this->getCluster(j).getRegion();
	Composition lo2 = clReg2.getOrigin();

	if (clReg1.isSimplex() and clReg2.isSimplex()) {
		// Look for H 1 and trap
		if ((lo1.isOnAxis(Species::H) and lo2.isOnAxis(Species::Trap)) or
			(lo2.isOnAxis(Species::H) and lo1.isOnAxis(Species::Trap))) {
			// Which one is which?
			auto trapId = lo1.isOnAxis(Species::Trap) ? i : j;
			auto hId = lo1.isOnAxis(Species::Trap) ? j : i;

			// Only single hydrogen can trap
			Composition loH = this->getCluster(hId).getRegion().getOrigin();
			if (loH[Species::H] == 1) {
				this->addTrapReaction(tag, {hId, trapId});
			}
		}
	}
}

inline ReactionCollection<LiReactionGenerator::NetworkType>
LiReactionGenerator::getReactionCollection() const
{
	ReactionCollection<NetworkType> ret(this->_clusterData.gridSize,
		this->_clusterData.numClusters, this->_enableReadRates,
		this->getProductionReactions(), this->getDissociationReactions(),
		this->getTrapReactions());
	return ret;
}
} // namespace detail

inline detail::LiReactionGenerator
LiReactionNetwork::getReactionGenerator() const noexcept
{
	return detail::LiReactionGenerator{*this};
}
} // namespace network
} // namespace core
} // namespace xolotl
