#pragma once

#include <algorithm>
#include <array>

#include <plsm/EnumIndexed.h>

#include <xolotl/core/Constants.h>
#include <xolotl/core/network/detail/ReactionUtility.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
Reaction<TNetwork, TDerived>::Reaction(ReactionDataRef reactionData,
	ClusterDataRef clusterData, IndexType reactionId) :
	_clusterData(clusterData),
	_reactionId(reactionId),
	_rate(reactionData.getRates(reactionId)),
	_widths(reactionData.getWidths(reactionId)),
	_coefs(reactionData.getCoefficients(reactionId))
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
Reaction<TNetwork, TDerived>::updateData(
	ReactionDataRef reactionData, ClusterDataRef clusterData)
{
	_clusterData = clusterData;
	_rate = reactionData.getRates(_reactionId);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
typename Reaction<TNetwork, TDerived>::AmountType
Reaction<TNetwork, TDerived>::computeOverlap(const ReflectedRegion& cl1RR,
	const ReflectedRegion& cl2RR, const ReflectedRegion& pr1RR,
	const ReflectedRegion& pr2RR)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	AmountType nOverlap = 1;
	for (auto i : speciesRangeNoI) {
		// The width is the subset of the tiles for which the
		// reaction is possible
		// For instance if we have X_1 + X_[3,5) ⇄ X_[5,7)
		// It is only possible for 4 within X_[3,5) and 5 within X_[5,7)
		// so the width is 1
		// More complicated with X_[3,5) + X_[5,7) ⇄ X_[9,11)
		// 3+6, 4+5, 4+6, width is 3

		for (auto m : makeIntervalRange(pr2RR[i()]))
			for (auto l : makeIntervalRange(cl2RR[i()])) {
				_widths(i()) += util::max(0.0,
					util::min(
						cl1RR[i()].end() - 1 + l, pr1RR[i()].end() - 1 + m) -
						util::max(
							cl1RR[i()].begin() + l, pr1RR[i()].begin() + m) +
						1.0);
			}

		nOverlap *= _widths(i());
	}

	//	if (nOverlap <= 0) {
	//		std::cout << "first reactant: ";
	//		for (auto i : speciesRangeNoI) {
	//			std::cout << cl1RR[i()].begin() << ", ";
	//		}
	//		std::cout << std::endl;
	//		for (auto i : speciesRangeNoI) {
	//			std::cout << cl1RR[i()].end() - 1 << ", ";
	//		}
	//		std::cout << std::endl << "second reactant: ";
	//		for (auto i : speciesRangeNoI) {
	//			std::cout << cl2RR[i()].begin() << ", ";
	//		}
	//		std::cout << std::endl;
	//		for (auto i : speciesRangeNoI) {
	//			std::cout << cl2RR[i()].end() - 1 << ", ";
	//		}
	//		std::cout << std::endl << "product: ";
	//		for (auto i : speciesRangeNoI) {
	//			std::cout << pr1RR[i()].begin() << ", ";
	//		}
	//		std::cout << std::endl;
	//		for (auto i : speciesRangeNoI) {
	//			std::cout << pr1RR[i()].end() - 1 << ", ";
	//		}
	//		std::cout << std::endl << "second product: ";
	//		for (auto i : speciesRangeNoI) {
	//			std::cout << pr2RR[i()].begin() << ", ";
	//		}
	//		std::cout << std::endl;
	//		for (auto i : speciesRangeNoI) {
	//			std::cout << pr2RR[i()].end() - 1 << ", ";
	//		}
	//		std::cout << std::endl;
	//		std::cout << "Overlap: " << nOverlap << std::endl;
	//		std::cout << "Widths: ";
	//		for (auto i : speciesRangeNoI) {
	//			std::cout << _widths(i()) << ", ";
	//		}
	//		std::cout << std::endl;
	//	}
	assert(nOverlap > 0);

	return nOverlap;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
ProductionReaction<TNetwork, TDerived>::ProductionReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, IndexType cluster0, IndexType cluster1,
	IndexType cluster2, IndexType cluster3) :
	Superclass(reactionData, clusterData, reactionId),
	_reactants({cluster0, cluster1}),
	_products({cluster2, cluster3})
{
	for (auto i : {0, 1}) {
		this->copyMomentIds(_reactants[i], _reactantMomentIds[i]);
		this->copyMomentIds(_products[i], _productMomentIds[i]);
	}

	this->initialize();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
ProductionReaction<TNetwork, TDerived>::ProductionReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	ProductionReaction(reactionData, clusterData, reactionId,
		clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2,
		clusterSet.cluster3)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeCoefficients()
{
	// static
	const auto dummyRegion = Region(Composition{});

	// Find the overlap for this reaction
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	const auto& cl1Reg =
		this->_clusterData.getCluster(_reactants[0]).getRegion();
	const auto& cl2Reg =
		this->_clusterData.getCluster(_reactants[1]).getRegion();
	const auto& pr1Reg = (_products[0] == invalidIndex) ?
		dummyRegion :
		this->_clusterData.getCluster(_products[0]).getRegion();
	const auto& pr2Reg = (_products[1] == invalidIndex) ?
		dummyRegion :
		this->_clusterData.getCluster(_products[1]).getRegion();
	const auto& cl1Disp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			cl1Reg);
	const auto& cl2Disp =
		detail::getReflectedDispersionForCoefs<NetworkType::Traits::numSpecies>(
			cl2Reg);

	// Initialize the reflected regions
	auto rRegions = detail::updateReflectedRegionsForCoefs<nMomentIds>(
		cl1Reg, cl2Reg, pr1Reg, pr2Reg);
	auto cl1RR = rRegions[0];
	auto cl2RR = rRegions[1];
	auto pr1RR = rRegions[2];
	auto pr2RR = rRegions[3];

	// If there is no product the overlap is 1
	double nOverlap = 1.0;
	// No product case
	if (_products[0] == invalidIndex && _products[1] == invalidIndex) {
		for (auto i : speciesRangeNoI) {
			this->_widths[i()] = 1.0;
		}
	}
	// General case
	else
		nOverlap = static_cast<double>(
			this->computeOverlap(cl1RR, cl2RR, pr1RR, pr2RR));

	this->_coefs(0, 0, 0, 0) = nOverlap;
	for (auto i : speciesRangeNoI) {
		// First order sum on the first reactant
		auto factor = nOverlap / this->_widths[i()];
		this->_coefs(i() + 1, 0, 0, 0) = factor *
			detail::computeFirstOrderSum(i(), cl1RR, cl2RR, pr1RR, pr2RR);

		this->_coefs(0, 0, 0, i() + 1) =
			this->_coefs(i() + 1, 0, 0, 0) / cl1Disp[i()];

		// First order sum on the second reactant
		this->_coefs(0, i() + 1, 0, 0) = factor *
			detail::computeFirstOrderSum(i(), cl2RR, cl1RR, pr1RR, pr2RR);

		this->_coefs(0, 0, 1, i() + 1) =
			this->_coefs(0, i() + 1, 0, 0) / cl2Disp[i()];

		// Loop on the potential products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}

			// Get the regions in the right order
			const auto& thisRR = (prodId == _products[0]) ? pr1RR : pr2RR;
			const auto& otherRR = (prodId == _products[0]) ? pr2RR : pr1RR;
			// Get the dispersion
			const auto& thisDispersion = (prodId == _products[0]) ?
				detail::getReflectedDispersionForCoefs<
					NetworkType::Traits::numSpecies>(pr1Reg) :
				detail::getReflectedDispersionForCoefs<
					NetworkType::Traits::numSpecies>(pr2Reg);

			// First order sum on the other product (p+2) because 0 and 1 are
			// used for reactants)
			this->_coefs(0, 0, p + 2, i() + 1) = factor *
				detail::computeFirstOrderSum(
					i(), thisRR, otherRR, cl2RR, cl1RR) /
				thisDispersion[i()];

			// Products first moments
			for (auto k : speciesRangeNoI) {
				// Second order sum
				if (k == i) {
					this->_coefs(i() + 1, 0, p + 2, k() + 1) = factor *
						detail::computeSecondOrderOffsetSum(
							i(), cl1RR, cl2RR, thisRR, otherRR) /
						thisDispersion[i()];

					this->_coefs(0, i() + 1, p + 2, k() + 1) = factor *
						detail::computeSecondOrderOffsetSum(
							i(), cl2RR, cl1RR, thisRR, otherRR) /
						thisDispersion[i()];
				}
				else {
					this->_coefs(i() + 1, 0, p + 2, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, 0, p + 2, k() + 1) / nOverlap;

					this->_coefs(0, i() + 1, p + 2, k() + 1) =
						this->_coefs(0, i() + 1, 0, 0) *
						this->_coefs(0, 0, p + 2, k() + 1) / nOverlap;
				}
			}
		}
	}

	for (auto i : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[i()];
		for (auto j : speciesRangeNoI) {
			// Second order sum
			if (i == j) {
				for (double m : makeIntervalRange(pr2RR[j()]))
					for (double l : makeIntervalRange(cl1RR[j()])) {
						this->_coefs(i() + 1, j() + 1, 0, 0) +=
							(l -
								static_cast<double>(
									cl1RR[j()].end() - 1 + cl1RR[j()].begin()) /
									2.0) *
							factor *
							util::firstOrderSum(
								util::max(pr1RR[j()].begin() + m - l,
									static_cast<double>(cl2RR[j()].begin())),
								util::min(pr1RR[j()].end() - 1 + m - l,
									static_cast<double>(cl2RR[j()].end() - 1)),
								static_cast<double>(
									cl2RR[j()].end() - 1 + cl2RR[j()].begin()) /
									2.0);
					}
			}
			else {
				this->_coefs(i() + 1, j() + 1, 0, 0) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(0, j() + 1, 0, 0) / nOverlap;
			}
		}
	}

	for (auto i : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[i()];

		// First reactant first moments
		for (auto k : speciesRangeNoI) {
			if (k == i) {
				this->_coefs(i() + 1, 0, 0, k() + 1) = factor *
					detail::computeSecondOrderSum(
						i(), cl1RR, cl2RR, pr1RR, pr2RR) /
					cl1Disp[i()];
			}
			else {
				this->_coefs(i() + 1, 0, 0, k() + 1) =
					this->_coefs(i() + 1, 0, 0, 0) *
					this->_coefs(k() + 1, 0, 0, 0) / (nOverlap * cl1Disp[k()]);
			}

			this->_coefs(0, i() + 1, 0, k() + 1) =
				this->_coefs(k() + 1, i() + 1, 0, 0) / cl1Disp[k()];
		}

		// Second reactant partial derivatives
		for (auto k : speciesRangeNoI) {
			if (k == i) {
				this->_coefs(0, i() + 1, 1, k() + 1) = factor *
					detail::computeSecondOrderSum(
						i(), cl2RR, cl1RR, pr1RR, pr2RR) /
					cl2Disp[i()];
			}
			else {
				this->_coefs(0, i() + 1, 1, k() + 1) =
					this->_coefs(0, i() + 1, 0, 0) *
					this->_coefs(0, k() + 1, 0, 0) / (nOverlap * cl2Disp[k()]);
			}

			this->_coefs(i() + 1, 0, 1, k() + 1) =
				this->_coefs(i() + 1, k() + 1, 0, 0) / cl2Disp[k()];
		}
	}

	// Now we loop over the 2 dimensions of the coefs to compute all
	// the possible sums over distances for the flux
	for (auto i : speciesRangeNoI) {
		auto factor = nOverlap / this->_widths[i()];
		for (auto j : speciesRangeNoI) {
			// Now we deal with the coefficients needed for the
			// first moments
			// Let's start with the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}

				// Get the regions in the right order
				const auto& thisRR = (prodId == _products[0]) ? pr1RR : pr2RR;
				const auto& otherRR = (prodId == _products[0]) ? pr2RR : pr1RR;
				// Get the dispersion
				const auto& thisDispersion = (prodId == _products[0]) ?
					detail::getReflectedDispersionForCoefs<
						NetworkType::Traits::numSpecies>(pr1Reg) :
					detail::getReflectedDispersionForCoefs<
						NetworkType::Traits::numSpecies>(pr2Reg);

				for (auto k : speciesRangeNoI) {
					// Third order sum
					if (i == j && j == k) {
						for (double m : makeIntervalRange(otherRR[i()]))
							for (double l : makeIntervalRange(cl1RR[i()])) {
								this->_coefs(
									i() + 1, j() + 1, p + 2, k() + 1) +=
									(l -
										static_cast<double>(cl1RR[i()].end() -
											1 + cl1RR[i()].begin()) /
											2.0) *
									factor *
									util::secondOrderOffsetSum(
										util::max(thisRR[i()].begin() + m - l,
											static_cast<double>(
												cl2RR[i()].begin())),
										util::min(thisRR[i()].end() - 1 + m - l,
											static_cast<double>(
												cl2RR[i()].end() - 1)),
										static_cast<double>(cl2RR[i()].end() -
											1 + cl2RR[i()].begin()) /
											2.0,
										static_cast<double>(thisRR[i()].end() -
											1 + thisRR[i()].begin()) /
											2.0,
										l - m);
							}
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) /=
							thisDispersion[k()];
					}
					else if (j == k) {
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) =
							this->_coefs(i() + 1, 0, 0, 0) *
							this->_coefs(0, j() + 1, p + 2, k() + 1) / nOverlap;
					}
					else if (i == k) {
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) =
							this->_coefs(0, j() + 1, 0, 0) *
							this->_coefs(i() + 1, 0, p + 2, k() + 1) / nOverlap;
					}
					else if (i == j) {
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) =
							this->_coefs(0, 0, p + 2, k() + 1) *
							this->_coefs(i() + 1, j() + 1, 0, 0) / nOverlap;
					}
					else {
						this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) =
							this->_coefs(i() + 1, 0, 0, 0) *
							this->_coefs(0, j() + 1, 0, 0) *
							this->_coefs(0, 0, p + 2, k() + 1) /
							(nOverlap * nOverlap);
					}
				}
			}

			// Let's take care of the first reactant first moments
			for (auto k : speciesRangeNoI) {
				// Third order sum
				if (i == j && j == k) {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) = factor *
						detail::computeThirdOrderSum(
							i(), cl2RR, cl1RR, pr1RR, pr2RR) /
						cl1Disp[i()];
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) = factor *
						detail::computeThirdOrderSum(
							i(), cl1RR, cl2RR, pr1RR, pr2RR) /
						cl2Disp[i()];
				}
				else if (i == k) {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) =
						this->_coefs(0, j() + 1, 0, 0) *
						this->_coefs(i() + 1, 0, 0, k() + 1) / nOverlap;
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) =
						this->_coefs(0, j() + 1, 0, 0) *
						this->_coefs(i() + 1, 0, 1, k() + 1) / nOverlap;
				}
				else if (j == k) {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, j() + 1, 0, k() + 1) / nOverlap;
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, j() + 1, 1, k() + 1) / nOverlap;
				}
				else if (i == j) {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) =
						this->_coefs(0, 0, 0, k() + 1) *
						this->_coefs(i() + 1, j() + 1, 0, 0) / nOverlap;
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) =
						this->_coefs(0, 0, 1, k() + 1) *
						this->_coefs(i() + 1, j() + 1, 0, 0) / nOverlap;
				}
				else {
					this->_coefs(i() + 1, j() + 1, 0, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, j() + 1, 0, 0) *
						this->_coefs(k() + 1, 0, 0, 0) /
						(nOverlap * nOverlap * cl1Disp[k()]);
					this->_coefs(i() + 1, j() + 1, 1, k() + 1) =
						this->_coefs(i() + 1, 0, 0, 0) *
						this->_coefs(0, j() + 1, 0, 0) *
						this->_coefs(0, k() + 1, 0, 0) /
						(nOverlap * nOverlap * cl2Disp[k()]);
				}
			}
		}
	}
}

template <typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() > 1),
	double>
getRateForProduction(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg,
	const double r0, const double r1, const double dc0, const double dc1)
{
	constexpr double pi = ::xolotl::core::pi;
	constexpr double rCore = ::xolotl::core::alloyCoreRadius;
	const double zs = 4.0 * pi * (r0 + r1 + rCore);
	using Species = typename TRegion::EnumIndex;
	bool cl0IsSphere = (pairCl0Reg.getOrigin().isOnAxis(Species::V) ||
			 pairCl0Reg.getOrigin().isOnAxis(Species::Void) ||
			 pairCl0Reg.getOrigin().isOnAxis(Species::I)),
		 cl1IsSphere = (pairCl1Reg.getOrigin().isOnAxis(Species::V) ||
			 pairCl1Reg.getOrigin().isOnAxis(Species::Void) ||
			 pairCl1Reg.getOrigin().isOnAxis(Species::I));

	// Simple case
	if (cl0IsSphere && cl1IsSphere) {
		return zs * (dc0 + dc1);
	}

	double p = 0.0, zl = 0.0;
	if (r0 < r1) {
		p = 1.0 / (1.0 + pow(r1 / (3.0 * (r0 + rCore)), 2.0));
		zl = 4.0 * pow(pi, 2.0) * r1 / log(1.0 + 8.0 * r1 / (r0 + rCore));
	}
	else {
		p = 1.0 / (1.0 + pow(r0 / (3.0 * (r1 + rCore)), 2.0));
		zl = 4.0 * pow(pi, 2.0) * r0 / log(1.0 + 8.0 * r0 / (r1 + rCore));
	}

	double k_plus = (dc0 + dc1) * (p * zs + (1.0 - p) * zl);
	double bias = 1.0;
	if (pairCl0Reg.getOrigin().isOnAxis(Species::I) ||
		pairCl1Reg.getOrigin().isOnAxis(Species::I)) {
		bias = 1.2;
	}

	return k_plus * bias;
}

template <typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() < 2),
	double>
getRateForProduction(const TRegion& pairCl0Reg, const TRegion& pairCl1Reg,
	const double r0, const double r1, const double dc0, const double dc1)
{
	constexpr double pi = ::xolotl::core::pi;

	double kPlus = 4.0 * pi * (r0 + r1) * (dc0 + dc1);

	return kPlus;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
ProductionReaction<TNetwork, TDerived>::computeRate(IndexType gridIndex)
{
	auto cl0 = this->_clusterData.getCluster(_reactants[0]);
	auto cl1 = this->_clusterData.getCluster(_reactants[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	return getRateForProduction(
		cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
	// Get the total number of elements in each cluster
	auto cl1 = this->_clusterData.getCluster(_reactants[0]);
	const auto& cl1Reg = cl1.getRegion();
	const bool cl1IsSimplex = cl1Reg.isSimplex();
	auto cl2 = this->_clusterData.getCluster(_reactants[1]);
	const auto& cl2Reg = cl2.getRegion();
	const bool cl2IsSimplex = cl2Reg.isSimplex();
	// Each reactant connects with all the reactants
	// Reactant 1 with reactant 1
	this->addConnectivity(_reactants[0], _reactants[0], connectivity);
	if (!cl1IsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_reactants[0], _reactantMomentIds[0][i()], connectivity);
			this->addConnectivity(
				_reactantMomentIds[0][i()], _reactants[0], connectivity);
			for (auto j : speciesRangeNoI) {
				this->addConnectivity(_reactantMomentIds[0][i()],
					_reactantMomentIds[0][j()], connectivity);
			}
		}
	}
	// Reactant 2 with reactant 1
	this->addConnectivity(_reactants[1], _reactants[0], connectivity);
	if (!cl1IsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_reactants[1], _reactantMomentIds[0][i()], connectivity);
		}
	}
	if (!cl2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_reactantMomentIds[1][i()], _reactants[0], connectivity);
		}
	}
	if (!cl1IsSimplex && !cl2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			for (auto j : speciesRangeNoI) {
				this->addConnectivity(_reactantMomentIds[1][i()],
					_reactantMomentIds[0][j()], connectivity);
			}
		}
	}
	// Reactant 1 with reactant 2
	this->addConnectivity(_reactants[0], _reactants[1], connectivity);
	if (!cl2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_reactants[0], _reactantMomentIds[1][i()], connectivity);
		}
	}
	if (!cl1IsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_reactantMomentIds[0][i()], _reactants[1], connectivity);
		}
	}
	if (!cl1IsSimplex && !cl2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			for (auto j : speciesRangeNoI) {
				this->addConnectivity(_reactantMomentIds[0][i()],
					_reactantMomentIds[1][j()], connectivity);
			}
		}
	}
	// Reactant 2 with reactant 2
	this->addConnectivity(_reactants[1], _reactants[1], connectivity);
	if (!cl2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_reactants[1], _reactantMomentIds[1][i()], connectivity);
			this->addConnectivity(
				_reactantMomentIds[1][i()], _reactants[1], connectivity);
			for (auto j : speciesRangeNoI) {
				this->addConnectivity(_reactantMomentIds[1][i()],
					_reactantMomentIds[1][j()], connectivity);
			}
		}
	}
	// Each product connects with all the reactants
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		auto prod = this->_clusterData.getCluster(prodId);
		const auto& prodReg = prod.getRegion();
		const bool prodIsSimplex = prodReg.isSimplex();

		// With reactant 1
		this->addConnectivity(prodId, _reactants[0], connectivity);
		if (!cl1IsSimplex) {
			for (auto i : speciesRangeNoI) {
				this->addConnectivity(
					prodId, _reactantMomentIds[0][i()], connectivity);
			}
		}
		if (!prodIsSimplex) {
			for (auto i : speciesRangeNoI) {
				this->addConnectivity(
					_productMomentIds[p][i()], _reactants[0], connectivity);
			}
		}
		if (!cl1IsSimplex && !prodIsSimplex) {
			for (auto i : speciesRangeNoI) {
				for (auto j : speciesRangeNoI) {
					this->addConnectivity(_productMomentIds[p][i()],
						_reactantMomentIds[0][j()], connectivity);
				}
			}
		}
		// With reactant 2
		this->addConnectivity(prodId, _reactants[1], connectivity);
		if (!cl2IsSimplex) {
			for (auto i : speciesRangeNoI) {
				this->addConnectivity(
					prodId, _reactantMomentIds[1][i()], connectivity);
			}
		}
		if (!prodIsSimplex) {
			for (auto i : speciesRangeNoI) {
				this->addConnectivity(
					_productMomentIds[p][i()], _reactants[1], connectivity);
			}
		}
		if (!cl2IsSimplex && !prodIsSimplex) {
			for (auto i : speciesRangeNoI) {
				for (auto j : speciesRangeNoI) {
					this->addConnectivity(_productMomentIds[p][i()],
						_reactantMomentIds[1][j()], connectivity);
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
	// Get the total number of elements in each cluster
	auto cl1 = this->_clusterData.getCluster(_reactants[0]);
	const auto& cl1Reg = cl1.getRegion();
	const bool cl1IsSimplex = cl1Reg.isSimplex();
	auto cl2 = this->_clusterData.getCluster(_reactants[1]);
	const auto& cl2Reg = cl2.getRegion();
	const bool cl2IsSimplex = cl2Reg.isSimplex();
	// Each reactant connects with all the reactants
	// Reactant 1 with reactant 1
	this->addConnectivity(_reactants[0], _reactants[0], connectivity);
	if (!cl1IsSimplex) {
		for (auto i : speciesRangeNoI) {
			for (auto j : speciesRangeNoI) {
				if (i() == j())
					this->addConnectivity(_reactantMomentIds[0][i()],
						_reactantMomentIds[0][j()], connectivity);
			}
		}
	}
	// Reactant 2 with reactant 1
	if (_reactants[1] == _reactants[0])
		this->addConnectivity(_reactants[1], _reactants[0], connectivity);
	if (!cl1IsSimplex && !cl2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[1][i()] == _reactantMomentIds[0][j()])
					this->addConnectivity(_reactantMomentIds[1][i()],
						_reactantMomentIds[0][j()], connectivity);
			}
		}
	}
	// Reactant 1 with reactant 2
	if (_reactants[1] == _reactants[0])
		this->addConnectivity(_reactants[0], _reactants[1], connectivity);
	if (!cl1IsSimplex && !cl2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			for (auto j : speciesRangeNoI) {
				if (_reactantMomentIds[0][i()] == _reactantMomentIds[1][j()])
					this->addConnectivity(_reactantMomentIds[0][i()],
						_reactantMomentIds[1][j()], connectivity);
			}
		}
	}
	// Reactant 2 with reactant 2
	this->addConnectivity(_reactants[1], _reactants[1], connectivity);
	if (!cl2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			for (auto j : speciesRangeNoI) {
				if (i() == j())
					this->addConnectivity(_reactantMomentIds[1][i()],
						_reactantMomentIds[1][j()], connectivity);
			}
		}
	}
	// Each product connects with all the reactants
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		auto prod = this->_clusterData.getCluster(prodId);
		const auto& prodReg = prod.getRegion();
		const bool prodIsSimplex = prodReg.isSimplex();

		// With reactant 1
		if (prodId == _reactants[0])
			this->addConnectivity(prodId, _reactants[0], connectivity);
		if (!cl1IsSimplex && !prodIsSimplex) {
			for (auto i : speciesRangeNoI) {
				for (auto j : speciesRangeNoI) {
					if (_productMomentIds[p][i()] == _reactantMomentIds[0][j()])
						this->addConnectivity(_productMomentIds[p][i()],
							_reactantMomentIds[0][j()], connectivity);
				}
			}
		}
		// With reactant 2
		if (prodId == _reactants[1])
			this->addConnectivity(prodId, _reactants[1], connectivity);
		if (!cl2IsSimplex && !prodIsSimplex) {
			for (auto i : speciesRangeNoI) {
				for (auto j : speciesRangeNoI) {
					if (_productMomentIds[p][i()] == _reactantMomentIds[1][j()])
						this->addConnectivity(_productMomentIds[p][i()],
							_reactantMomentIds[1][j()], connectivity);
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeFlux(
	ConcentrationsView concentrations, FluxesView fluxes, IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Compute the total number of elements in each cluster
	auto cl1 = this->_clusterData.getCluster(_reactants[0]);
	const auto& cl1Reg = cl1.getRegion();
	AmountType volCl1 = cl1Reg.volume();
	auto cl2 = this->_clusterData.getCluster(_reactants[1]);
	const auto& cl2Reg = cl2.getRegion();
	AmountType volCl2 = cl2Reg.volume();

	// Compute the flux for the 0th order moments
	double f = this->_coefs(0, 0, 0, 0) * concentrations[_reactants[0]] *
		concentrations[_reactants[1]];
	for (auto i : speciesRangeNoI) {
		f += this->_coefs(i() + 1, 0, 0, 0) *
			concentrations[_reactantMomentIds[0][i()]] *
			concentrations[_reactants[1]];
	}
	for (auto j : speciesRangeNoI) {
		f += this->_coefs(0, j() + 1, 0, 0) * concentrations[_reactants[0]] *
			concentrations[_reactantMomentIds[1][j()]];
	}
	for (auto i : speciesRangeNoI) {
		for (auto j : speciesRangeNoI) {
			f += this->_coefs(i() + 1, j() + 1, 0, 0) *
				concentrations[_reactantMomentIds[0][i()]] *
				concentrations[_reactantMomentIds[1][j()]];
		}
	}
	f *= this->_rate(gridIndex);

	Kokkos::atomic_sub(&fluxes[_reactants[0]], f / (double)volCl1);
	Kokkos::atomic_sub(&fluxes[_reactants[1]], f / (double)volCl2);
	for (auto prodId : _products) {
		if (prodId == invalidIndex) {
			continue;
		}

		auto prod = this->_clusterData.getCluster(prodId);
		const auto& prodReg = prod.getRegion();
		AmountType volProd = prodReg.volume();
		Kokkos::atomic_add(&fluxes[prodId], f / (double)volProd);
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		// First for the first reactant
		if (volCl1 > 1) {
			f = this->_coefs(0, 0, 0, k() + 1) * concentrations[_reactants[0]] *
				concentrations[_reactants[1]];
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 0, k() + 1) *
					concentrations[_reactantMomentIds[0][i()]] *
					concentrations[_reactants[1]];
			}
			for (auto j : speciesRangeNoI) {
				f += this->_coefs(0, j() + 1, 0, k() + 1) *
					concentrations[_reactants[0]] *
					concentrations[_reactantMomentIds[1][j()]];
			}
			for (auto i : speciesRangeNoI) {
				for (auto j : speciesRangeNoI) {
					f += this->_coefs(i() + 1, j() + 1, 0, k() + 1) *
						concentrations[_reactantMomentIds[0][i()]] *
						concentrations[_reactantMomentIds[1][j()]];
				}
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_sub(
				&fluxes[_reactantMomentIds[0][k()]], f / (double)volCl1);
		}

		// For the second reactant
		if (volCl2 > 1) {
			f = this->_coefs(0, 0, 1, k() + 1) * concentrations[_reactants[0]] *
				concentrations[_reactants[1]];
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 1, k() + 1) *
					concentrations[_reactantMomentIds[0][i()]] *
					concentrations[_reactants[1]];
			}
			for (auto j : speciesRangeNoI) {
				f += this->_coefs(0, j() + 1, 1, k() + 1) *
					concentrations[_reactants[0]] *
					concentrations[_reactantMomentIds[1][j()]];
			}
			for (auto i : speciesRangeNoI) {
				for (auto j : speciesRangeNoI) {
					f += this->_coefs(i() + 1, j() + 1, 1, k() + 1) *
						concentrations[_reactantMomentIds[0][i()]] *
						concentrations[_reactantMomentIds[1][j()]];
				}
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_sub(
				&fluxes[_reactantMomentIds[1][k()]], f / (double)volCl2);
		}

		// For the products
		for (auto p : {0, 1}) {
			auto prodId = _products[p];
			if (prodId == invalidIndex) {
				continue;
			}

			auto prod = this->_clusterData.getCluster(prodId);
			const auto& prodReg = prod.getRegion();
			AmountType volProd = prodReg.volume();

			if (volProd > 1) {
				f = this->_coefs(0, 0, p + 2, k() + 1) *
					concentrations[_reactants[0]] *
					concentrations[_reactants[1]];
				for (auto i : speciesRangeNoI) {
					f += this->_coefs(i() + 1, 0, p + 2, k() + 1) *
						concentrations[_reactantMomentIds[0][i()]] *
						concentrations[_reactants[1]];
				}
				for (auto j : speciesRangeNoI) {
					f += this->_coefs(0, j() + 1, p + 2, k() + 1) *
						concentrations[_reactants[0]] *
						concentrations[_reactantMomentIds[1][j()]];
				}
				for (auto i : speciesRangeNoI) {
					for (auto j : speciesRangeNoI) {
						f += this->_coefs(i() + 1, j() + 1, p + 2, k() + 1) *
							concentrations[_reactantMomentIds[0][i()]] *
							concentrations[_reactantMomentIds[1][j()]];
					}
				}
				f *= this->_rate(gridIndex);
				Kokkos::atomic_add(
					&fluxes[_productMomentIds[p][k()]], f / (double)volProd);
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computePartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
	int nProd = 0;
	for (auto prodId : _products) {
		if (prodId != invalidIndex) {
			++nProd;
		}
	}

	// Compute the total number of elements in each cluster
	auto cl1 = this->_clusterData.getCluster(_reactants[0]);
	const auto& cl1Reg = cl1.getRegion();
	AmountType volCl1 = cl1Reg.volume();
	auto cl2 = this->_clusterData.getCluster(_reactants[1]);
	const auto& cl2Reg = cl2.getRegion();
	AmountType volCl2 = cl2Reg.volume();

	// Compute the partials for the 0th order moments
	// Compute the values (d / dL_0^A)
	double temp = this->_coefs(0, 0, 0, 0) * concentrations[_reactants[1]];
	if (volCl2 > 1) {
		for (auto i : speciesRangeNoI) {
			temp += this->_coefs(0, i() + 1, 0, 0) *
				concentrations[_reactantMomentIds[1][i()]];
		}
	}
	// First for the first reactant
	Kokkos::atomic_sub(&values(connectivity(_reactants[0], _reactants[0])),
		this->_rate(gridIndex) * temp / (double)volCl1);
	// Second reactant
	Kokkos::atomic_sub(&values(connectivity(_reactants[1], _reactants[0])),
		this->_rate(gridIndex) * temp / (double)volCl2);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		auto prod = this->_clusterData.getCluster(prodId);
		const auto& prodReg = prod.getRegion();
		AmountType volProd = prodReg.volume();
		Kokkos::atomic_add(&values(connectivity(prodId, _reactants[0])),
			this->_rate(gridIndex) * temp / (double)volProd);
	}

	// Compute the values (d / dL_0^B)
	temp = this->_coefs(0, 0, 0, 0) * concentrations[_reactants[0]];
	if (volCl1 > 1) {
		for (auto i : speciesRangeNoI) {
			temp += this->_coefs(i() + 1, 0, 0, 0) *
				concentrations[_reactantMomentIds[0][i()]];
		}
	}
	// First for the first reactant
	Kokkos::atomic_sub(&values(connectivity(_reactants[0], _reactants[1])),
		this->_rate(gridIndex) * temp / (double)volCl1);
	// Second reactant
	Kokkos::atomic_sub(&values(connectivity(_reactants[1], _reactants[1])),
		this->_rate(gridIndex) * temp / (double)volCl2);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}
		auto prod = this->_clusterData.getCluster(prodId);
		const auto& prodReg = prod.getRegion();
		AmountType volProd = prodReg.volume();
		Kokkos::atomic_add(&values(connectivity(prodId, _reactants[1])),
			this->_rate(gridIndex) * temp / (double)volProd);
	}

	// (d / dL_1^A)
	if (volCl1 > 1) {
		for (auto i : speciesRangeNoI) {
			temp =
				this->_coefs(i() + 1, 0, 0, 0) * concentrations[_reactants[1]];
			if (volCl2 > 1) {
				for (auto j : speciesRangeNoI) {
					temp += this->_coefs(i() + 1, j() + 1, 0, 0) *
						concentrations[_reactantMomentIds[1][j()]];
				}
			}
			// First reactant
			Kokkos::atomic_sub(&values(connectivity(
								   _reactants[0], _reactantMomentIds[0][i()])),
				this->_rate(gridIndex) * temp / (double)volCl1);
			// second reactant
			Kokkos::atomic_sub(&values(connectivity(
								   _reactants[1], _reactantMomentIds[0][i()])),
				this->_rate(gridIndex) * temp / (double)volCl2);
			// For the products
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				auto prod = this->_clusterData.getCluster(prodId);
				const auto& prodReg = prod.getRegion();
				AmountType volProd = prodReg.volume();
				Kokkos::atomic_add(
					&values(connectivity(prodId, _reactantMomentIds[0][i()])),
					this->_rate(gridIndex) * temp / (double)volProd);
			}
		}
	}

	// (d / dL_1^B)
	if (volCl2 > 1) {
		for (auto i : speciesRangeNoI) {
			temp =
				this->_coefs(0, i() + 1, 0, 0) * concentrations[_reactants[0]];
			if (volCl1 > 1) {
				for (auto j : speciesRangeNoI) {
					temp += this->_coefs(j() + 1, i() + 1, 0, 0) *
						concentrations[_reactantMomentIds[0][j()]];
				}
			}
			Kokkos::atomic_sub(&values(connectivity(
								   _reactants[0], _reactantMomentIds[1][i()])),
				this->_rate(gridIndex) * temp / (double)volCl1);
			Kokkos::atomic_sub(&values(connectivity(
								   _reactants[1], _reactantMomentIds[1][i()])),
				this->_rate(gridIndex) * temp / (double)volCl2);
			for (auto p : {0, 1}) {
				auto prodId = _products[p];
				if (prodId == invalidIndex) {
					continue;
				}
				auto prod = this->_clusterData.getCluster(prodId);
				const auto& prodReg = prod.getRegion();
				AmountType volProd = prodReg.volume();
				Kokkos::atomic_add(
					&values(connectivity(prodId, _reactantMomentIds[1][i()])),
					this->_rate(gridIndex) * temp / (double)volProd);
			}
		}
	}

	// Take care of the first moments
	if (volCl1 > 1) {
		for (auto k : speciesRangeNoI) {
			// First for the first reactant
			// (d / dL_0^A)
			temp =
				this->_coefs(0, 0, 0, k() + 1) * concentrations[_reactants[1]];
			if (volCl2 > 1) {
				for (auto j : speciesRangeNoI) {
					temp += this->_coefs(0, j() + 1, 0, k() + 1) *
						concentrations[_reactantMomentIds[1][j()]];
				}
			}
			Kokkos::atomic_sub(&values(connectivity(
								   _reactantMomentIds[0][k()], _reactants[0])),
				this->_rate(gridIndex) * temp / (double)volCl1);
			// (d / dL_1^A)
			for (auto i : speciesRangeNoI) {
				temp = this->_coefs(i() + 1, 0, 0, k() + 1) *
					concentrations[_reactants[1]];
				if (volCl2 > 1) {
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(i() + 1, j() + 1, 0, k() + 1) *
							concentrations[_reactantMomentIds[1][j()]];
					}
				}
				Kokkos::atomic_sub(
					&values(connectivity(_reactantMomentIds[0][k()],
						_reactantMomentIds[0][i()])),
					this->_rate(gridIndex) * temp / (double)volCl1);
			}
			// (d / dL_0^B)
			temp =
				this->_coefs(0, 0, 0, k() + 1) * concentrations[_reactants[0]];
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(j() + 1, 0, 0, k() + 1) *
					concentrations[_reactantMomentIds[0][j()]];
			}
			Kokkos::atomic_sub(&values(connectivity(
								   _reactantMomentIds[0][k()], _reactants[1])),
				this->_rate(gridIndex) * temp / (double)volCl1);
			// (d / dL_1^B)
			if (volCl2 > 1) {
				for (auto i : speciesRangeNoI) {
					temp = this->_coefs(0, i() + 1, 0, k() + 1) *
						concentrations[_reactants[0]];
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(j() + 1, i() + 1, 0, k() + 1) *
							concentrations[_reactantMomentIds[0][j()]];
					}
					Kokkos::atomic_sub(
						&values(connectivity(_reactantMomentIds[0][k()],
							_reactantMomentIds[1][i()])),
						this->_rate(gridIndex) * temp / (double)volCl1);
				}
			}
		}
	}

	// Take care of the first moments
	if (volCl2 > 1) {
		for (auto k : speciesRangeNoI) {
			// First for the second reactant
			// (d / dL_0^A)
			temp =
				this->_coefs(0, 0, 1, k() + 1) * concentrations[_reactants[1]];
			for (auto j : speciesRangeNoI) {
				temp += this->_coefs(0, j() + 1, 1, k() + 1) *
					concentrations[_reactantMomentIds[1][j()]];
			}
			Kokkos::atomic_sub(&values(connectivity(
								   _reactantMomentIds[1][k()], _reactants[0])),
				this->_rate(gridIndex) * temp / (double)volCl2);
			// (d / dL_1^A)
			if (volCl1 > 1) {
				for (auto i : speciesRangeNoI) {
					temp = this->_coefs(i() + 1, 0, 1, k() + 1) *
						concentrations[_reactants[1]];
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(i() + 1, j() + 1, 1, k() + 1) *
							concentrations[_reactantMomentIds[1][j()]];
					}
					Kokkos::atomic_sub(
						&values(connectivity(_reactantMomentIds[1][k()],
							_reactantMomentIds[0][i()])),
						this->_rate(gridIndex) * temp / (double)volCl2);
				}
			}
			// (d / dL_0^B)
			temp =
				this->_coefs(0, 0, 1, k() + 1) * concentrations[_reactants[0]];
			if (volCl1 > 1) {
				for (auto j : speciesRangeNoI) {
					temp += this->_coefs(j() + 1, 0, 1, k() + 1) *
						concentrations[_reactantMomentIds[0][j()]];
				}
			}
			Kokkos::atomic_sub(&values(connectivity(
								   _reactantMomentIds[1][k()], _reactants[1])),
				this->_rate(gridIndex) * temp / (double)volCl2);
			// (d / dL_1^B)
			for (auto i : speciesRangeNoI) {
				temp = this->_coefs(0, i() + 1, 1, k() + 1) *
					concentrations[_reactants[0]];
				if (volCl1 > 1) {
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(j() + 1, i() + 1, 1, k() + 1) *
							concentrations[_reactantMomentIds[0][j()]];
					}
				}
				Kokkos::atomic_sub(
					&values(connectivity(_reactantMomentIds[1][k()],
						_reactantMomentIds[1][i()])),
					this->_rate(gridIndex) * temp / (double)volCl2);
			}
		}
	}

	// Loop on the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}

		auto prod = this->_clusterData.getCluster(prodId);
		const auto& prodReg = prod.getRegion();
		AmountType volProd = prodReg.volume();

		// Take care of the first moments
		if (volProd > 1) {
			for (auto k : speciesRangeNoI) {
				// (d / dL_0^A)
				temp = this->_coefs(0, 0, p + 2, k() + 1) *
					concentrations[_reactants[1]];
				if (volCl2 > 1) {
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(0, j() + 1, p + 2, k() + 1) *
							concentrations[_reactantMomentIds[1][j()]];
					}
				}
				Kokkos::atomic_add(
					&values(
						connectivity(_productMomentIds[p][k()], _reactants[0])),
					this->_rate(gridIndex) * temp / (double)volProd);
				// (d / dL_1^A)
				if (volCl1 > 1) {
					for (auto i : speciesRangeNoI) {
						temp = this->_coefs(i() + 1, 0, p + 2, k() + 1) *
							concentrations[_reactants[1]];
						if (volCl2 > 1) {
							for (auto j : speciesRangeNoI) {
								temp += this->_coefs(
											i() + 1, j() + 1, p + 2, k() + 1) *
									concentrations[_reactantMomentIds[1][j()]];
							}
						}
						Kokkos::atomic_add(
							&values(connectivity(_productMomentIds[p][k()],
								_reactantMomentIds[0][i()])),
							this->_rate(gridIndex) * temp / (double)volProd);
					}
				}
				// (d / dL_0^B)
				temp = this->_coefs(0, 0, p + 2, k() + 1) *
					concentrations[_reactants[0]];
				if (volCl1 > 1) {
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(j() + 1, 0, p + 2, k() + 1) *
							concentrations[_reactantMomentIds[0][j()]];
					}
				}
				Kokkos::atomic_add(
					&values(
						connectivity(_productMomentIds[p][k()], _reactants[1])),
					this->_rate(gridIndex) * temp / (double)volProd);
				// (d / dL_1^B)
				if (volCl2 > 1) {
					for (auto i : speciesRangeNoI) {
						temp = this->_coefs(0, i() + 1, p + 2, k() + 1) *
							concentrations[_reactants[0]];
						if (volCl1 > 1) {
							for (auto j : speciesRangeNoI) {
								temp += this->_coefs(
											j() + 1, i() + 1, p + 2, k() + 1) *
									concentrations[_reactantMomentIds[0][j()]];
							}
						}
						Kokkos::atomic_add(
							&values(connectivity(_productMomentIds[p][k()],
								_reactantMomentIds[1][i()])),
							this->_rate(gridIndex) * temp / (double)volProd);
					}
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ProductionReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
	ConcentrationsView concentrations, Kokkos::View<double*> values,
	Connectivity connectivity, IndexType gridIndex)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();
	int nProd = 0;
	for (auto prodId : _products) {
		if (prodId != invalidIndex) {
			++nProd;
		}
	}

	// Compute the total number of elements in each cluster
	auto cl1 = this->_clusterData.getCluster(_reactants[0]);
	const auto& cl1Reg = cl1.getRegion();
	AmountType volCl1 = cl1Reg.volume();
	auto cl2 = this->_clusterData.getCluster(_reactants[1]);
	const auto& cl2Reg = cl2.getRegion();
	AmountType volCl2 = cl2Reg.volume();

	// Compute the partials for the 0th order moments
	// Compute the values (d / dL_0^A)
	double temp = this->_coefs(0, 0, 0, 0) * concentrations[_reactants[1]];
	if (volCl2 > 1) {
		for (auto i : speciesRangeNoI) {
			temp += this->_coefs(0, i() + 1, 0, 0) *
				concentrations[_reactantMomentIds[1][i()]];
		}
	}
	// First for the first reactant
	Kokkos::atomic_sub(&values(connectivity(_reactants[0], _reactants[0])),
		this->_rate(gridIndex) * temp / (double)volCl1);
	// Second reactant
	if (_reactants[1] == _reactants[0])
		Kokkos::atomic_sub(&values(connectivity(_reactants[1], _reactants[0])),
			this->_rate(gridIndex) * temp / (double)volCl2);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex || prodId != _reactants[0]) {
			continue;
		}
		auto prod = this->_clusterData.getCluster(prodId);
		const auto& prodReg = prod.getRegion();
		AmountType volProd = prodReg.volume();
		Kokkos::atomic_add(&values(connectivity(prodId, _reactants[0])),
			this->_rate(gridIndex) * temp / (double)volProd);
	}

	// Compute the values (d / dL_0^B)
	temp = this->_coefs(0, 0, 0, 0) * concentrations[_reactants[0]];
	if (volCl1 > 1) {
		for (auto i : speciesRangeNoI) {
			temp += this->_coefs(i() + 1, 0, 0, 0) *
				concentrations[_reactantMomentIds[0][i()]];
		}
	}
	// First for the first reactant
	if (_reactants[1] == _reactants[0])
		Kokkos::atomic_sub(&values(connectivity(_reactants[0], _reactants[1])),
			this->_rate(gridIndex) * temp / (double)volCl1);
	// Second reactant
	Kokkos::atomic_sub(&values(connectivity(_reactants[1], _reactants[1])),
		this->_rate(gridIndex) * temp / (double)volCl2);
	// For the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex || prodId != _reactants[1]) {
			continue;
		}
		auto prod = this->_clusterData.getCluster(prodId);
		const auto& prodReg = prod.getRegion();
		AmountType volProd = prodReg.volume();
		Kokkos::atomic_add(&values(connectivity(prodId, _reactants[1])),
			this->_rate(gridIndex) * temp / (double)volProd);
	}

	// Take care of the first moments
	if (volCl1 > 1) {
		for (auto k : speciesRangeNoI) {
			// First for the first reactant
			// (d / dL_1^A)
			for (auto i : speciesRangeNoI) {
				temp = this->_coefs(i() + 1, 0, 0, k() + 1) *
					concentrations[_reactants[1]];
				if (volCl2 > 1) {
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(i() + 1, j() + 1, 0, k() + 1) *
							concentrations[_reactantMomentIds[1][j()]];
					}
				}
				if (k() == i())
					Kokkos::atomic_sub(
						&values(connectivity(_reactantMomentIds[0][k()],
							_reactantMomentIds[0][i()])),
						this->_rate(gridIndex) * temp / (double)volCl1);
			}
			// (d / dL_1^B)
			if (volCl2 > 1) {
				for (auto i : speciesRangeNoI) {
					temp = this->_coefs(0, i() + 1, 0, k() + 1) *
						concentrations[_reactants[0]];
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(j() + 1, i() + 1, 0, k() + 1) *
							concentrations[_reactantMomentIds[0][j()]];
					}
					if (_reactantMomentIds[0][k()] ==
						_reactantMomentIds[1][i()])
						Kokkos::atomic_sub(
							&values(connectivity(_reactantMomentIds[0][k()],
								_reactantMomentIds[1][i()])),
							this->_rate(gridIndex) * temp / (double)volCl1);
				}
			}
		}
	}

	// Take care of the first moments
	if (volCl2 > 1) {
		for (auto k : speciesRangeNoI) {
			// First for the second reactant
			// (d / dL_1^A)
			if (volCl1 > 1) {
				for (auto i : speciesRangeNoI) {
					temp = this->_coefs(i() + 1, 0, 1, k() + 1) *
						concentrations[_reactants[1]];
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(i() + 1, j() + 1, 1, k() + 1) *
							concentrations[_reactantMomentIds[1][j()]];
					}
					if (_reactantMomentIds[1][k()] ==
						_reactantMomentIds[0][i()])
						Kokkos::atomic_sub(
							&values(connectivity(_reactantMomentIds[1][k()],
								_reactantMomentIds[0][i()])),
							this->_rate(gridIndex) * temp / (double)volCl2);
				}
			}
			// (d / dL_1^B)
			for (auto i : speciesRangeNoI) {
				temp = this->_coefs(0, i() + 1, 1, k() + 1) *
					concentrations[_reactants[0]];
				if (volCl1 > 1) {
					for (auto j : speciesRangeNoI) {
						temp += this->_coefs(j() + 1, i() + 1, 1, k() + 1) *
							concentrations[_reactantMomentIds[0][j()]];
					}
				}
				if (k() == i())
					Kokkos::atomic_sub(
						&values(connectivity(_reactantMomentIds[1][k()],
							_reactantMomentIds[1][i()])),
						this->_rate(gridIndex) * temp / (double)volCl2);
			}
		}
	}

	// Loop on the products
	for (auto p : {0, 1}) {
		auto prodId = _products[p];
		if (prodId == invalidIndex) {
			continue;
		}

		auto prod = this->_clusterData.getCluster(prodId);
		const auto& prodReg = prod.getRegion();
		AmountType volProd = prodReg.volume();

		// Take care of the first moments
		if (volProd > 1) {
			for (auto k : speciesRangeNoI) {
				// (d / dL_1^A)
				if (volCl1 > 1) {
					for (auto i : speciesRangeNoI) {
						temp = this->_coefs(i() + 1, 0, p + 2, k() + 1) *
							concentrations[_reactants[1]];
						if (volCl2 > 1) {
							for (auto j : speciesRangeNoI) {
								temp += this->_coefs(
											i() + 1, j() + 1, p + 2, k() + 1) *
									concentrations[_reactantMomentIds[1][j()]];
							}
						}
						if (_productMomentIds[p][k()] ==
							_reactantMomentIds[0][i()])
							Kokkos::atomic_add(
								&values(connectivity(_productMomentIds[p][k()],
									_reactantMomentIds[0][i()])),
								this->_rate(gridIndex) * temp /
									(double)volProd);
					}
				}
				// (d / dL_1^B)
				if (volCl2 > 1) {
					for (auto i : speciesRangeNoI) {
						temp = this->_coefs(0, i() + 1, p + 2, k() + 1) *
							concentrations[_reactants[0]];
						if (volCl1 > 1) {
							for (auto j : speciesRangeNoI) {
								temp += this->_coefs(
											j() + 1, i() + 1, p + 2, k() + 1) *
									concentrations[_reactantMomentIds[0][j()]];
							}
						}
						if (_productMomentIds[p][k()] ==
							_reactantMomentIds[1][i()])
							Kokkos::atomic_add(
								&values(connectivity(_productMomentIds[p][k()],
									_reactantMomentIds[1][i()])),
								this->_rate(gridIndex) * temp /
									(double)volProd);
					}
				}
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
ProductionReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	// Check if our cluster is on the left side of this reaction
	if (clusterId == _reactants[0]) {
		return this->_rate(gridIndex) * concentrations[_reactants[1]] *
			this->_coefs(0, 0, 0, 0);
	}
	if (clusterId == _reactants[1]) {
		return this->_rate(gridIndex) * concentrations[_reactants[0]] *
			this->_coefs(0, 0, 0, 0);
	}

	// This cluster is not part of the reaction
	return 0.0;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
DissociationReaction<TNetwork, TDerived>::DissociationReaction(
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
DissociationReaction<TNetwork, TDerived>::DissociationReaction(
	ReactionDataRef reactionData, ClusterDataRef clusterData,
	IndexType reactionId, const detail::ClusterSet& clusterSet) :
	DissociationReaction(reactionData, clusterData, reactionId,
		clusterSet.cluster0, clusterSet.cluster1, clusterSet.cluster2)
{
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeCoefficients()
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
DissociationReaction<TNetwork, TDerived>::computeRate(IndexType gridIndex)
{
	double omega = this->_clusterData.getAtomicVolume();
	double T = this->_clusterData.temperature(gridIndex);

	// TODO: computeProductionRate should use products and not reactants
	auto cl0 = this->_clusterData.getCluster(_products[0]);
	auto cl1 = this->_clusterData.getCluster(_products[1]);

	double r0 = cl0.getReactionRadius();
	double r1 = cl1.getReactionRadius();

	double dc0 = cl0.getDiffusionCoefficient(gridIndex);
	double dc1 = cl1.getDiffusionCoefficient(gridIndex);

	double kPlus = getRateForProduction(
		cl0.getRegion(), cl1.getRegion(), r0, r1, dc0, dc1);
	double E_b = this->asDerived()->computeBindingEnergy();

	constexpr double k_B = ::xolotl::core::kBoltzmann;

	double kMinus = (1.0 / omega) * kPlus * std::exp(-E_b / (k_B * T));

	return kMinus;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Get the total number of elements in each cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	const bool clIsSimplex = clReg.isSimplex();
	auto prod1 = this->_clusterData.getCluster(_products[0]);
	const auto& prod1Reg = prod1.getRegion();
	const bool prod1IsSimplex = prod1Reg.isSimplex();
	auto prod2 = this->_clusterData.getCluster(_products[1]);
	const auto& prod2Reg = prod2.getRegion();
	const bool prod2IsSimplex = prod2Reg.isSimplex();

	// The reactant connects with the reactant
	this->addConnectivity(_reactant, _reactant, connectivity);
	if (!clIsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_reactant, _reactantMomentIds[i()], connectivity);
			this->addConnectivity(
				_reactantMomentIds[i()], _reactant, connectivity);
			for (auto j : speciesRangeNoI) {
				this->addConnectivity(_reactantMomentIds[i()],
					_reactantMomentIds[j()], connectivity);
			}
		}
	}
	// Each product connects with  the reactant
	// Product 1 with reactant
	this->addConnectivity(_products[0], _reactant, connectivity);
	if (!clIsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_products[0], _reactantMomentIds[i()], connectivity);
		}
	}
	if (!prod1IsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_productMomentIds[0][i()], _reactant, connectivity);
		}
	}
	if (!clIsSimplex && !prod1IsSimplex) {
		for (auto i : speciesRangeNoI) {
			for (auto j : speciesRangeNoI) {
				this->addConnectivity(_productMomentIds[0][i()],
					_reactantMomentIds[j()], connectivity);
			}
		}
	}
	// Product 2 with reactant
	this->addConnectivity(_products[1], _reactant, connectivity);
	if (!clIsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_products[1], _reactantMomentIds[i()], connectivity);
		}
	}
	if (!prod2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			this->addConnectivity(
				_productMomentIds[1][i()], _reactant, connectivity);
		}
	}
	if (!clIsSimplex && !prod2IsSimplex) {
		for (auto i : speciesRangeNoI) {
			for (auto j : speciesRangeNoI) {
				this->addConnectivity(_productMomentIds[1][i()],
					_reactantMomentIds[j()], connectivity);
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeReducedConnectivity(
	const Connectivity& connectivity)
{
	constexpr auto speciesRangeNoI = NetworkType::getSpeciesRangeNoI();

	// Get the total number of elements in each cluster
	auto cl = this->_clusterData.getCluster(_reactant);
	const auto& clReg = cl.getRegion();
	const bool clIsSimplex = clReg.isSimplex();
	auto prod1 = this->_clusterData.getCluster(_products[0]);
	const auto& prod1Reg = prod1.getRegion();
	const bool prod1IsSimplex = prod1Reg.isSimplex();
	auto prod2 = this->_clusterData.getCluster(_products[1]);
	const auto& prod2Reg = prod2.getRegion();
	const bool prod2IsSimplex = prod2Reg.isSimplex();

	// The reactant connects with the reactant
	this->addConnectivity(_reactant, _reactant, connectivity);
	if (!clIsSimplex) {
		for (auto i : speciesRangeNoI) {
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
	if (!clIsSimplex && !prod1IsSimplex) {
		for (auto i : speciesRangeNoI) {
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
	if (!clIsSimplex && !prod2IsSimplex) {
		for (auto i : speciesRangeNoI) {
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
DissociationReaction<TNetwork, TDerived>::computeFlux(
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

	// Compute the flux for the 0th order moments
	double f = this->_coefs(0, 0, 0, 0) * concentrations[_reactant];
	for (auto i : speciesRangeNoI) {
		f += this->_coefs(i() + 1, 0, 0, 0) *
			concentrations[_reactantMomentIds[i()]];
	}
	f *= this->_rate(gridIndex);
	Kokkos::atomic_sub(&fluxes[_reactant], f / (double)volCl);
	Kokkos::atomic_add(&fluxes[_products[0]], f / (double)volProd1);
	Kokkos::atomic_add(&fluxes[_products[1]], f / (double)volProd2);

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		// First for the reactant
		if (volCl > 1) {
			f = this->_coefs(0, 0, 0, k() + 1) * concentrations[_reactant];
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 0, k() + 1) *
					concentrations[_reactantMomentIds[i()]];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_sub(
				&fluxes[_reactantMomentIds[k()]], f / (double)volCl);
		}

		// Now the first product
		if (volProd1 > 1) {
			f = this->_coefs(0, 0, 1, k() + 1) * concentrations[_reactant];
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 1, k() + 1) *
					concentrations[_reactantMomentIds[i()]];
			}
			f *= this->_rate(gridIndex);
			Kokkos::atomic_add(
				&fluxes[_productMomentIds[0][k()]], f / (double)volProd1);
		}

		// Finally the second product
		if (volProd2 > 1) {
			f = this->_coefs(0, 0, 2, k() + 1) * concentrations[_reactant];
			for (auto i : speciesRangeNoI) {
				f += this->_coefs(i() + 1, 0, 2, k() + 1) *
					concentrations[_reactantMomentIds[i()]];
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
DissociationReaction<TNetwork, TDerived>::computePartialDerivatives(
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
	if (volProd1 > 1) {
		for (auto i : speciesRangeNoI) {
			Kokkos::atomic_sub(
				&values(connectivity(_reactant, _reactantMomentIds[i()])),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}
	// For the first product
	df = this->_rate(gridIndex) / (double)volProd1;
	Kokkos::atomic_add(&values(connectivity(_products[0], _reactant)),
		df * this->_coefs(0, 0, 0, 0));

	if (volProd1 > 1) {
		for (auto i : speciesRangeNoI) {
			Kokkos::atomic_add(
				&values(connectivity(_products[0], _reactantMomentIds[i()])),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}
	// For the second product
	df = this->_rate(gridIndex) / (double)volProd2;
	Kokkos::atomic_add(&values(connectivity(_products[1], _reactant)),
		df * this->_coefs(0, 0, 0, 0));

	if (volProd1 > 1) {
		for (auto i : speciesRangeNoI) {
			Kokkos::atomic_add(
				&values(connectivity(_products[1], _reactantMomentIds[i()])),
				df * this->_coefs(i() + 1, 0, 0, 0));
		}
	}

	// Take care of the first moments
	for (auto k : speciesRangeNoI) {
		if (volCl > 1) {
			// First for the reactant
			df = this->_rate(gridIndex) / (double)volCl;
			// Compute the values
			Kokkos::atomic_sub(
				&values(connectivity(_reactantMomentIds[k()], _reactant)),
				df * this->_coefs(0, 0, 0, k() + 1));
			for (auto i : speciesRangeNoI) {
				Kokkos::atomic_sub(&values(connectivity(_reactantMomentIds[k()],
									   _reactantMomentIds[i()])),
					df * this->_coefs(i() + 1, 0, 0, k() + 1));
			}
		}
		// For the first product
		if (volProd1 > 1) {
			df = this->_rate(gridIndex) / (double)volProd1;
			Kokkos::atomic_add(
				&values(connectivity(_productMomentIds[0][k()], _reactant)),
				df * this->_coefs(0, 0, 1, k() + 1));
			for (auto i : speciesRangeNoI) {
				Kokkos::atomic_add(
					&values(connectivity(
						_productMomentIds[0][k()], _reactantMomentIds[i()])),
					df * this->_coefs(i() + 1, 0, 1, k() + 1));
			}
		}
		// For the second product
		if (volProd2 > 1) {
			df = this->_rate(gridIndex) / (double)volProd2;
			Kokkos::atomic_add(
				&values(connectivity(_productMomentIds[1][k()], _reactant)),
				df * this->_coefs(0, 0, 2, k() + 1));
			for (auto i : speciesRangeNoI) {
				Kokkos::atomic_add(
					&values(connectivity(
						_productMomentIds[1][k()], _reactantMomentIds[i()])),
					df * this->_coefs(i() + 1, 0, 2, k() + 1));
			}
		}
	}
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
DissociationReaction<TNetwork, TDerived>::computeReducedPartialDerivatives(
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
		if (volCl > 1) {
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
		if (volProd1 > 1) {
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
		if (volProd2 > 1) {
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
DissociationReaction<TNetwork, TDerived>::computeLeftSideRate(
	ConcentrationsView concentrations, IndexType clusterId, IndexType gridIndex)
{
	// Check if our cluster is on the left side of this reaction
	if (clusterId == _reactant) {
		return this->_rate(gridIndex) * this->_coefs(0, 0, 0, 0);
	}

	// This cluster is not part of the reaction
	return 0.0;
}
} // namespace network
} // namespace core
} // namespace xolotl
