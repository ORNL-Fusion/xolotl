#pragma once

#include <xolotl/core/network/detail/ReactionUtility.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
ReSolutionReaction<TNetwork, TDerived>::ReSolutionReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1,
	IndexType cluster2) :
	Superclass(reactionData, clusterData, reactionId),
	_reactant(cluster0),
	_products({cluster1, cluster2})
{
	this->copyMomentIds(_reactant, _reactantMomentIds);
	for (auto i : {0, 1}) {
		this->copyMomentIds(_products[i], _productMomentIds[i]);
	}

	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
ReSolutionReaction<TNetwork, TDerived>::ReSolutionReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	ReSolutionReaction(reactionData, clusterData, reactionId,
		clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReSolutionReaction<TNetwork, TDerived>::computeCoefficients()
{
	// static
	const auto dummyRegion = Region(Composition{});

	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	auto clReg = this->_clusterData.getCluster(_reactant).getRegion();
	auto prod1Reg = this->_clusterData.getCluster(_products[0]).getRegion();
	auto prod2Reg = this->_clusterData.getCluster(_products[1]).getRegion();
	const auto& clDisp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			clReg);
	const auto& prod1Disp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			prod1Reg);
	const auto& prod2Disp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			prod2Reg);
	auto cl2Reg = dummyRegion;

	// Initialize the reflected regions
	auto rRegions = detail::updateReflectedRegionsForCoefs<nMomentIds>(
		prod1Reg, prod2Reg, clReg, cl2Reg);
	auto clRR = rRegions[2];
	auto cl2RR = rRegions[3];
	auto pr1RR = rRegions[0];
	auto pr2RR = rRegions[1];

	auto nOverlap =
		static_cast<double>(this->computeOverlap(pr1RR, pr2RR, clRR, cl2RR));

	// The first coefficient is simply the overlap because it is the sum over 1
	this->_coefs(0, 0, 0, 0) = nOverlap;
	for (auto i : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[i()];
		// First order sum
		this->_coefs(i() + 1, 0, 0, 0) = factor *
			detail::computeFirstOrderSum(i(), clRR, cl2RR, pr2RR, pr1RR);
	}

	// First moments
	for (auto k : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[k()];
		// Reactant
		this->_coefs(0, 0, 0, k() + 1) =
			this->_coefs(k() + 1, 0, 0, 0) / clDisp[k()];

		// First product
		this->_coefs(0, 0, 1, k() + 1) = factor *
			detail::computeFirstOrderSum(k(), pr1RR, pr2RR, clRR, cl2RR) /
			prod1Disp[k()];

		// Second product
		this->_coefs(0, 0, 2, k() + 1) = factor *
			detail::computeFirstOrderSum(k(), pr2RR, pr1RR, clRR, cl2RR) /
			prod2Disp[k()];
	}

	// Now we loop over the 1 dimension of the coefs to compute all the
	// possible sums over distances for the flux
	for (auto i : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[i()];
		// Now we deal with the coefficients needed for the partial derivatives
		// Starting with the reactant
		for (auto k : speciesRangeNoI) {
			// Second order sum
			if (k == i) {
				this->_coefs(i() + 1, 0, 0, k() + 1) = factor *
					detail::computeSecondOrderSum(
						i(), clRR, cl2RR, pr2RR, pr1RR) /
					clDisp[k()];
			}
			else {
				this->_coefs(i() + 1, 0, 0, k() + 1) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(k() + 1, 0, 0, 0) / (nOverlap * clDisp[k()]);
			}
		}

		// First moments for the first product
		for (auto k : speciesRangeNoI) {
			// Second order sum
			if (k == i) {
				this->_coefs(i() + 1, 0, 1, k() + 1) = factor *
					detail::computeSecondOrderOffsetSum(
						i(), clRR, cl2RR, pr1RR, pr2RR) /
					prod1Disp[k()];
			}
			else {
				this->_coefs(i() + 1, 0, 1, k() + 1) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(0, 0, 1, k() + 1) / nOverlap;
			}
		}

		// First moments for the second product
		for (auto k : speciesRangeNoI) {
			// Second order sum
			if (k == i) {
				this->_coefs(i() + 1, 0, 2, k() + 1) = factor *
					detail::computeSecondOrderOffsetSum(
						i(), clRR, cl2RR, pr2RR, pr1RR) /
					prod2Disp[k()];
			}
			else {
				this->_coefs(i() + 1, 0, 2, k() + 1) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(0, 0, 2, k() + 1) / nOverlap;
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
ReSolutionReaction<TNetwork, TDerived>::computeRate(IndexType gridIndex)
{
	// Get Zeta
	auto zeta = this->_clusterData.zeta(0);
	// Set the fit variables depending on the electronic stopping power (zeta)
	double y0 = 0.0, a1 = 0.0, b1 = 0.0, b2 = 0.0, c = 0.0;
	if (zeta > 0.87) {
		y0 = 11.0851, a1 = 1.5052, b1 = 0.0362, b2 = 0.0203, c = 3.4123;
	}
	else if (zeta > 0.82) {
		y0 = 10.6297, a1 = 1.3479, b1 = 0.0438, b2 = 0.0241, c = 4.2214;
	}
	else if (zeta > 0.77) {
		y0 = 10.1521, a1 = 1.1986, b1 = 0.0546, b2 = 0.0299, c = 5.4612;
	}
	else if (zeta > 0.71) {
		y0 = 9.1816, a1 = 0.949, b1 = 0.0703, b2 = 0.0371, c = 7.982;
	}
	else if (zeta > 0.67) {
		y0 = 8.6745, a1 = 0.8401, b1 = 0.0792, b2 = 0.0407, c = 9.6585;
	}
	else if (zeta > 0.62) {
		y0 = 7.6984, a1 = 0.6721, b1 = 0.1028, b2 = 0.0526, c = 14.272;
	}
	else if (zeta > 0.57) {
		y0 = 6.3925, a1 = 0.5025, b1 = 0.1411, b2 = 0.0727, c = 23.1967;
	}
	else if (zeta > 0.52) {
		y0 = 4.6175, a1 = 0.3433, b1 = 0.2284, b2 = 0.1276, c = 45.6624;
	}
	else {
		y0 = 2.3061, a1 = 0.2786, b1 = 1.1008, b2 = 1.605, c = 150.6689;
	}
	// Compute the fraction rate
	auto cl = this->_clusterData.getCluster(_reactant);
	auto radius = cl.getReactionRadius();
	double fractionRate = (a1 * exp(-b1 * radius) +
							  (y0 - a1) / (1.0 + c * pow(radius, 2.0)) *
								  exp(-b2 * pow(radius, 2.0))) *
		1.0e-4;
	// Resolution rate depends on fission rate
	double resolutionRate = 1.0e8 * this->_clusterData.fissionRate(0);

	return fractionRate * resolutionRate;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReSolutionReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Get the total number of elements in each cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	auto prod1 = this->_clusterData.getCluster(_products[0]);
	const auto& prod1Reg = prod1.getRegion();
	auto prod2 = this->_clusterData.getCluster(_products[1]);
	const auto& prod2Reg = prod2.getRegion();

	// The reactant connects with the reactant
	this->addConnectivity(_reactant, _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			this->addConnectivity(
				_reactant, _reactantMomentIds[i()], connectivity);
			this->addConnectivity(
				_reactantMomentIds[i()], _reactant, connectivity);
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[j()] != invalidIndex) {
					this->addConnectivity(_reactantMomentIds[i()],
						_reactantMomentIds[j()], connectivity);
				}
			}
		}
	}
	// Each product connects with  the reactant
	// Product 1 with reactant
	this->addConnectivity(_products[0], _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			this->addConnectivity(
				_products[0], _reactantMomentIds[i()], connectivity);
		}
	}
	for (auto i : speciesRangeNoI) {
		if (_productMomentIds[0][i()] != invalidIndex) {
			this->addConnectivity(
				_productMomentIds[0][i()], _reactant, connectivity);
		}
	}
	for (auto i : speciesRangeNoI) {
		if (_productMomentIds[0][i()] != invalidIndex) {
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[j()] != invalidIndex) {
					this->addConnectivity(_productMomentIds[0][i()],
						_reactantMomentIds[j()], connectivity);
				}
			}
		}
	}
	// Product 2 with reactant
	this->addConnectivity(_products[1], _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			this->addConnectivity(
				_products[1], _reactantMomentIds[i()], connectivity);
		}
	}
	for (auto i : speciesRangeNoI) {
		if (_productMomentIds[1][i()] != invalidIndex) {
			this->addConnectivity(
				_productMomentIds[1][i()], _reactant, connectivity);
		}
	}
	for (auto i : speciesRangeNoI) {
		if (_productMomentIds[1][i()] != invalidIndex) {
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[j()] != invalidIndex) {
					this->addConnectivity(_productMomentIds[1][i()],
						_reactantMomentIds[j()], connectivity);
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReSolutionReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Get the total number of elements in each cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	auto prod1 = this->_clusterData.getCluster(_products[0]);
	const auto& prod1Reg = prod1.getRegion();
	auto prod2 = this->_clusterData.getCluster(_products[1]);
	const auto& prod2Reg = prod2.getRegion();

	// The reactant connects with the reactant
	this->addConnectivity(_reactant, _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			for (auto j : speciesRangeNoI) {
				if (i() == j())
					this->addConnectivity(_reactantMomentIds[i()],
						_reactantMomentIds[j()], connectivity);
			}
		}
	}
	// Each product connects with  the reactant
	// Product 1 with reactant
	if (_products[0] == _reactant)
		this->addConnectivity(_products[0], _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_productMomentIds[0][i()] != invalidIndex) {
			for (auto j : speciesRangeNoI) {
				if (_productMomentIds[0][i()] == _reactantMomentIds[j()])
					this->addConnectivity(_productMomentIds[0][i()],
						_reactantMomentIds[j()], connectivity);
			}
		}
	}
	// Product 2 with reactant
	if (_products[1] == _reactant)
		this->addConnectivity(_products[1], _reactant, connectivity);
	for (auto i : speciesRangeNoI) {
		if (_productMomentIds[1][i()] != invalidIndex) {
			for (auto j : speciesRangeNoI) {
				if (_productMomentIds[1][i()] == _reactantMomentIds[j()])
					this->addConnectivity(_productMomentIds[1][i()],
						_reactantMomentIds[j()], connectivity);
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReSolutionReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the total number of elements in each cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	AmountType volCl = clReg.volume();
	auto prod1 = this->_clusterData.getCluster(_products[0]);
	const auto& prod1Reg = prod1.getRegion();
	AmountType volProd1 = prod1Reg.volume();
	auto prod2 = this->_clusterData.getCluster(_products[1]);
	const auto& prod2Reg = prod2.getRegion();
	AmountType volProd2 = prod2Reg.volume();

	// Initialize the concentrations that will be used in the loops
	auto cR = concentrations[_reactant];
	Kokkos::Array<double, nMomentIds> cmR;
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] == invalidIndex) {
			cmR[i()] = 0.0;
		}
		else
			cmR[i()] = concentrations[_reactantMomentIds[i()]];
	}

	// Compute the flux for the 0th order moments
	double f = this->_coefs(0, 0, 0, 0) * cR;
	for (auto i : speciesRangeNoI) {
		f += this->_coefs(i() + 1, 0, 0, 0) * cmR[i()];
	}
	f *= this->_rate(gridIndex);
	Kokkos::atomic_sub(&fluxes[_reactant], f / (double)volCl);
	Kokkos::atomic_add(&fluxes[_products[0]], f / (double)volProd1);
	Kokkos::atomic_add(&fluxes[_products[1]], f / (double)volProd2);

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		// First for the reactant
		if (_reactantMomentIds[k()] != invalidIndex) {
			f = this->_coefs(0, 0, 0, k() + 1) * cR;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 0, k() + 1) * cmR[i()];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_sub(
				&fluxes[_reactantMomentIds[k()]], f / (double)volCl);
		}

		// Now the first product
		if (_productMomentIds[0][k()] != invalidIndex) {
			f = this->_coefs(0, 0, 1, k() + 1) * cR;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 1, k() + 1) * cmR[i()];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_add(
				&fluxes[_productMomentIds[0][k()]], f / (double)volProd1);
		}

		// Finally the second product
		if (_productMomentIds[1][k()] != invalidIndex) {
			f = this->_coefs(0, 0, 2, k() + 1) * cR;
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 2, k() + 1) * cmR[i()];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_add(
				&fluxes[_productMomentIds[1][k()]], f / (double)volProd2);
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReSolutionReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	using AmountType = typename NetworkType::AmountType;
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the total number of elements in each cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	AmountType volCl = clReg.volume();
	auto prod1 = this->_clusterData.getCluster(_products[0]);
	const auto& prod1Reg = prod1.getRegion();
	AmountType volProd1 = prod1Reg.volume();
	auto prod2 = this->_clusterData.getCluster(_products[1]);
	const auto& prod2Reg = prod2.getRegion();
	AmountType volProd2 = prod2Reg.volume();

	// Compute the partials for the 0th order moments
	// First for the reactant
	double df = this->_rate(gridIndex) / (double)volCl;
	// Compute the values
	Kokkos::atomic_sub(&values(connectivity(_reactant, _reactant)),
		df * this->_coefs(0, 0, 0, 0));
	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_sub(
				&values(connectivity(_reactant, _reactantMomentIds[i()])),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}
	// For the first product
	df = this->_rate(gridIndex) / (double)volProd1;
	Kokkos::atomic_add(&values(connectivity(_products[0], _reactant)),
		df * this->_coefs(0, 0, 0, 0));

	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_add(
				&values(connectivity(_products[0], _reactantMomentIds[i()])),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}
	// For the second product
	df = this->_rate(gridIndex) / (double)volProd2;
	Kokkos::atomic_add(&values(connectivity(_products[1], _reactant)),
		df * this->_coefs(0, 0, 0, 0));

	for (auto i : speciesRangeNoI) {
		if (_reactantMomentIds[i()] != invalidIndex) {
			Kokkos::atomic_add(
				&values(connectivity(_products[1], _reactantMomentIds[i()])),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[k()] != invalidIndex) {
			// First for the reactant
			df = this->_rate(gridIndex) / (double)volCl;
			// Compute the values
			Kokkos::atomic_sub(
				&values(connectivity(_reactantMomentIds[k()], _reactant)),
				df * this->_coefs(0, 0, 0, k() + 1));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[i()] != invalidIndex) {
					Kokkos::atomic_sub(
						&values(connectivity(
							_reactantMomentIds[k()], _reactantMomentIds[i()])),
						df * this->_coefs(i() + 1, 0, 0, k() + 1));
				}
			}
		}
		// For the first product
		if (_productMomentIds[0][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / (double)volProd1;
			Kokkos::atomic_add(
				&values(connectivity(_productMomentIds[0][k()], _reactant)),
				df * this->_coefs(0, 0, 1, k() + 1));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[i()] != invalidIndex) {
					Kokkos::atomic_add(
						&values(connectivity(_productMomentIds[0][k()],
							_reactantMomentIds[i()])),
						df * this->_coefs(i() + 1, 0, 1, k() + 1));
				}
			}
		}
		// For the second product
		if (_productMomentIds[1][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / (double)volProd2;
			Kokkos::atomic_add(
				&values(connectivity(_productMomentIds[1][k()], _reactant)),
				df * this->_coefs(0, 0, 2, k() + 1));
			for (auto i : speciesRangeNoI) {
				if (_reactantMomentIds[i()] != invalidIndex) {
					Kokkos::atomic_add(
						&values(connectivity(_productMomentIds[1][k()],
							_reactantMomentIds[i()])),
						df * this->_coefs(i() + 1, 0, 2, k() + 1));
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReSolutionReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	using AmountType = typename NetworkType::AmountType;
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the total number of elements in each cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	AmountType volCl = clReg.volume();
	auto prod1 = this->_clusterData.getCluster(_products[0]);
	const auto& prod1Reg = prod1.getRegion();
	AmountType volProd1 = prod1Reg.volume();
	auto prod2 = this->_clusterData.getCluster(_products[1]);
	const auto& prod2Reg = prod2.getRegion();
	AmountType volProd2 = prod2Reg.volume();

	// Compute the partials for the 0th order moments
	// First for the reactant
	double df = this->_rate(gridIndex) / (double)volCl;
	// Compute the values
	Kokkos::atomic_sub(&values(connectivity(_reactant, _reactant)),
		df * this->_coefs(0, 0, 0, 0));
	// For the first product
	df = this->_rate(gridIndex) / (double)volProd1;
	if (_products[0] == _reactant)
		Kokkos::atomic_add(&values(connectivity(_products[0], _reactant)),
			df * this->_coefs(0, 0, 0, 0));

	// For the second product
	df = this->_rate(gridIndex) / (double)volProd2;
	if (_products[1] == _reactant)
		Kokkos::atomic_add(&values(connectivity(_products[1], _reactant)),
			df * this->_coefs(0, 0, 0, 0));

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (_reactantMomentIds[k()] != invalidIndex) {
			// First for the reactant
			df = this->_rate(gridIndex) / (double)volCl;
			// Compute the values
			for (auto i : speciesRangeNoI) {
				if (k() == i())
					Kokkos::atomic_sub(
						&values(connectivity(
							_reactantMomentIds[k()], _reactantMomentIds[i()])),
						df * this->_coefs(i() + 1, 0, 0, k() + 1));
			}
		}
		// For the first product
		if (_productMomentIds[0][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / (double)volProd1;
			for (auto i : speciesRangeNoI) {
				if (_productMomentIds[0][k()] == _reactantMomentIds[i()])
					Kokkos::atomic_add(
						&values(connectivity(_productMomentIds[0][k()],
							_reactantMomentIds[i()])),
						df * this->_coefs(i() + 1, 0, 1, k() + 1));
			}
		}
		// For the second product
		if (_productMomentIds[0][k()] != invalidIndex) {
			df = this->_rate(gridIndex) / (double)volProd2;
			for (auto i : speciesRangeNoI) {
				if (_productMomentIds[1][k()] == _reactantMomentIds[i()])
					Kokkos::atomic_add(
						&values(connectivity(_productMomentIds[1][k()],
							_reactantMomentIds[i()])),
						df * this->_coefs(i() + 1, 0, 2, k() + 1));
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
ReSolutionReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	// This type of reaction doesn't count
	return 0.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
