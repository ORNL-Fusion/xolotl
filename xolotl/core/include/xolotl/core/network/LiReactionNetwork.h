#pragma once

#include <xolotl/core/network/LiReaction.h>
#include <xolotl/core/network/LiTraits.h>
#include <xolotl/core/network/ReactionNetwork.h>
#include <xolotl/util/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class LiReactionGenerator;
}

class LiReactionNetwork : public ReactionNetwork<LiReactionNetwork>
{
	friend class ReactionNetwork<LiReactionNetwork>;

public:
	using Superclass = ReactionNetwork<LiReactionNetwork>;
	using Subpaving = typename Superclass::Subpaving;
	using Composition = typename Superclass::Composition;
	using Species = typename Superclass::Species;
	using AmountType = typename Superclass::AmountType;
	using IndexType = typename Superclass::IndexType;
	using ConcentrationsView = typename Superclass::ConcentrationsView;
	using FluxesView = typename Superclass::FluxesView;

	using Superclass::Superclass;

	IndexType
	checkLargestClusterId();

	void
	setTotalDepth(double depth) override
	{
		totalDepth = depth;
	}

	KOKKOS_INLINE_FUNCTION
	void setConnectivity(Connectivity);

	double
	computeTritiumRelease(
		double thetaH, double thetaT, double TSphereConc, double temperature);

	double
	getTritiumFlux(
		double thetaH, double thetaT, double TSphereConc, double temperature);

	double
	computeHydrogenPurge(double temperature)
	{
		double volumeFraction = 0.001; // 0.1 % H2
		double pressure = 1.0e5; // 1 atm

		return 2.0 * volumeFraction * pressure * nAv * 1.0e-27 /
			(gasConstant * temperature); // convert to nm-3
	}

	double
	getRT(double temperature)
	{
		return gasConstant * temperature * 1.0e18; // convert to nm2
	}

	double
	computeTritiumPurge(
		double thetaH, double thetaT, double hPurge, double temperature)
	{
		// Get all the terms
		auto aTerm = computeATritiumPurge(thetaH, thetaT, temperature);
		auto bTerm = computeBTritiumPurge(thetaH, thetaT, hPurge, temperature);
		auto cTerm = computeCTritiumPurge(thetaH, thetaT, hPurge, temperature);

		return (-bTerm + sqrt(bTerm * bTerm - 4.0 * aTerm * cTerm)) /
			(2.0 * aTerm);
	}

	//! Heat of adsorption, kg nm^2 s^-2 mol^-1
	double
	getQads(double theta)
	{
		if (theta > 0.1)
			return 1.0e18 * (20000.0 * (1.0 + pow(-log10(theta), 3.5)));
		return 1.0e18 * (-69999.108866 * log10(theta) - 29854.02);
	}

private:
	double
	checkLatticeParameter(double latticeParameter);

	double
	computeAtomicVolume(double latticeParameter)
	{
		// TODO
		return 0.5 * latticeParameter * latticeParameter * latticeParameter;
	}

	double
	checkImpurityRadius(double impurityRadius);

	detail::LiReactionGenerator
	getReactionGenerator() const noexcept;

	void
	readClusters(const std::string filename)
	{
		return;
	}

	void
	readReactions(double temperature, const std::string filename)
	{
		return;
	}

	void
	initializeExtraDOFs(const options::IOptions& options);

	void
	computeFluxesPreProcess(ConcentrationsView concentrations,
		FluxesView fluxes, IndexType gridIndex, double surfaceDepth,
		double spacing);

	void
	computePartialsPreProcess(ConcentrationsView concentrations,
		Kokkos::View<double*> values, IndexType gridIndex, double surfaceDepth,
		double spacing);

	double
	getQadsTheta(double theta)
	{
		if (theta > 0.1)
			return 1.0e18 *
				(-3.5 * 20000.0 * (pow(-log10(theta), 2.5)) /
					(theta * log(10)));
		return 1.0e18 * (-69999.108866 / (theta * log(10)));
	}

	double
	computeATritiumPurge(double thetaH, double thetaT, double temperature)
	{
		double RT = getRT(temperature);

		double t2AdsorptionRate =
			(condensationCoef * adjSites / sqrt(8.0e-3 * ::xolotl::core::pi)) *
			sqrt(RT / t2MolWeight) * exp(-2.0 * adsorptionE / RT);

		return (2.0 * purgeFlow) / sphereSurface +
			(t2AdsorptionRate * (1.0 - thetaH - thetaT) *
				(1.0 - thetaH - thetaT) * gamma_T2_T / 2.0);
	}

	double
	computeATritiumPurgeTheta(double thetaH, double thetaT, double temperature)
	{
		double RT = getRT(temperature);

		double t2AdsorptionRate =
			(condensationCoef * adjSites / sqrt(8.0e-3 * ::xolotl::core::pi)) *
			sqrt(RT / t2MolWeight) * exp(-2.0 * adsorptionE / RT);

		return -gamma_T2_T * t2AdsorptionRate * (1.0 - thetaH - thetaT);
	}

	double
	computeBTritiumPurge(
		double thetaH, double thetaT, double hPurge, double temperature)
	{
		double RT = getRT(temperature);

		double htAdsorptionRate =
			(condensationCoef * adjSites / sqrt(8.0e-3 * ::xolotl::core::pi)) *
			sqrt(RT / htMolWeight) * exp(-2.0 * adsorptionE / RT);

		double desorptionE = getQads(thetaH + thetaT) + adsorptionE;
		double desorptionRate =
			(surfaceSites * adjSites * RT * exp(-2.0 * desorptionE / RT)) /
			(2.0 * nAv * planckConstant);

		return hPurge * (2.0 * purgeFlow) / sphereSurface +
			(gamma_HT_T * htAdsorptionRate * (1.0 - thetaH - thetaT) *
				(1.0 - thetaH - thetaT) * hPurge) -
			(desorptionRate * (gamma_HT_T * thetaH + gamma_T2_T * thetaT) *
				thetaT);
	}

	double
	computeBTritiumPurgeThetaH(
		double thetaH, double thetaT, double hPurge, double temperature)
	{
		double RT = getRT(temperature);

		double htAdsorptionRate =
			(condensationCoef * adjSites / sqrt(8.0e-3 * ::xolotl::core::pi)) *
			sqrt(RT / htMolWeight) * exp(-2.0 * adsorptionE / RT);

		double desorptionE = getQads(thetaH + thetaT) + adsorptionE;
		double desorptionRate =
			(surfaceSites * adjSites * RT * exp(-2.0 * desorptionE / RT)) /
			(2.0 * nAv * planckConstant);
		double desorptionRateTheta = -surfaceSites * adjSites *
			getQadsTheta(thetaH + thetaT) * exp(-2.0 * desorptionE / RT) /
			(nAv * planckConstant);

		return -2.0 * gamma_HT_T * (1.0 - thetaH - thetaT) * htAdsorptionRate *
			hPurge -
			desorptionRateTheta * (gamma_HT_T * thetaH + gamma_T2_T * thetaT) *
			thetaT -
			desorptionRate * gamma_HT_T * thetaT;
	}

	double
	computeBTritiumPurgeThetaT(
		double thetaH, double thetaT, double hPurge, double temperature)
	{
		double RT = getRT(temperature);

		double htAdsorptionRate =
			(condensationCoef * adjSites / sqrt(8.0e-3 * ::xolotl::core::pi)) *
			sqrt(RT / htMolWeight) * exp(-2.0 * adsorptionE / RT);

		double desorptionE = getQads(thetaH + thetaT) + adsorptionE;
		double desorptionRate =
			(surfaceSites * adjSites * RT * exp(-2.0 * desorptionE / RT)) /
			(2.0 * nAv * planckConstant);
		double desorptionRateTheta = -surfaceSites * adjSites *
			getQadsTheta(thetaH + thetaT) * exp(-2.0 * desorptionE / RT) /
			(nAv * planckConstant);

		return -2.0 * gamma_HT_T * (1.0 - thetaH - thetaT) * htAdsorptionRate *
			hPurge -
			desorptionRateTheta * (gamma_HT_T * thetaH + gamma_T2_T * thetaT) *
			thetaT -
			desorptionRate * gamma_T2_T * thetaT -
			desorptionRate * (gamma_HT_T * thetaH + gamma_T2_T * thetaT);
	}

	double
	computeCTritiumPurge(
		double thetaH, double thetaT, double hPurge, double temperature)
	{
		double RT = getRT(temperature);

		double desorptionE = getQads(thetaH + thetaT) + adsorptionE;
		double desorptionRate =
			(surfaceSites * adjSites * RT * exp(-2.0 * desorptionE / RT)) /
			(2.0 * nAv * planckConstant);

		return -hPurge * desorptionRate *
			(gamma_HT_T * thetaH + gamma_T2_T * thetaT) * thetaT;
	}

	double
	computeCTritiumPurgeThetaH(
		double thetaH, double thetaT, double hPurge, double temperature)
	{
		double RT = getRT(temperature);

		double desorptionE = getQads(thetaH + thetaT) + adsorptionE;
		double desorptionRate =
			(surfaceSites * adjSites * RT * exp(-2.0 * desorptionE / RT)) /
			(2.0 * nAv * planckConstant);
		double desorptionRateTheta = -surfaceSites * adjSites *
			getQadsTheta(thetaH + thetaT) * exp(-2.0 * desorptionE / RT) /
			(nAv * planckConstant);

		return -desorptionRateTheta *
			(gamma_HT_T * thetaH + gamma_T2_T * thetaT) * thetaT * hPurge -
			desorptionRate * gamma_HT_T * thetaT * hPurge;
	}

	double
	computeCTritiumPurgeThetaT(
		double thetaH, double thetaT, double hPurge, double temperature)
	{
		double RT = getRT(temperature);

		double desorptionE = getQads(thetaH + thetaT) + adsorptionE;
		double desorptionRate =
			(surfaceSites * adjSites * RT * exp(-2.0 * desorptionE / RT)) /
			(2.0 * nAv * planckConstant);
		double desorptionRateTheta = -surfaceSites * adjSites *
			getQadsTheta(thetaH + thetaT) * exp(-2.0 * desorptionE / RT) /
			(nAv * planckConstant);

		return -desorptionRateTheta *
			(gamma_HT_T * thetaH + gamma_T2_T * thetaT) * thetaT * hPurge -
			desorptionRate * gamma_HT_T * thetaT * hPurge -
			desorptionRate * (gamma_HT_T * thetaH + gamma_T2_T * thetaT) *
			hPurge;
	}

	double
	computeTritiumPurgeThetaH(
		double thetaH, double thetaT, double hPurge, double temperature)
	{
		// Get all the terms
		auto aTerm = computeATritiumPurge(thetaH, thetaT, temperature);
		auto bTerm = computeBTritiumPurge(thetaH, thetaT, hPurge, temperature);
		auto cTerm = computeCTritiumPurge(thetaH, thetaT, hPurge, temperature);
		auto aTermTheta =
			computeATritiumPurgeTheta(thetaH, thetaT, temperature);
		auto bTermTheta =
			computeBTritiumPurgeThetaH(thetaH, thetaT, hPurge, temperature);
		auto cTermTheta =
			computeCTritiumPurgeThetaH(thetaH, thetaT, hPurge, temperature);

		double g = 2.0 * aTerm;
		double gPrime = 2.0 * aTermTheta;
		double h = bTerm * bTerm - 4.0 * aTerm * cTerm;
		double hPrime = 2.0 * bTermTheta * bTerm - 4.0 * aTermTheta * cTerm -
			4.0 * aTerm * cTermTheta;
		double f = -bTerm + sqrt(h);
		double fPrime = -bTermTheta + hPrime / (2.0 * sqrt(h));

		return (fPrime * g - f * gPrime) / (g * g);
	}

	double
	computeTritiumPurgeThetaT(
		double thetaH, double thetaT, double hPurge, double temperature)
	{
		// Get all the terms
		auto aTerm = computeATritiumPurge(thetaH, thetaT, temperature);
		auto bTerm = computeBTritiumPurge(thetaH, thetaT, hPurge, temperature);
		auto cTerm = computeCTritiumPurge(thetaH, thetaT, hPurge, temperature);
		auto aTermTheta =
			computeATritiumPurgeTheta(thetaH, thetaT, temperature);
		auto bTermTheta =
			computeBTritiumPurgeThetaT(thetaH, thetaT, hPurge, temperature);
		auto cTermTheta =
			computeCTritiumPurgeThetaT(thetaH, thetaT, hPurge, temperature);

		double g = 2.0 * aTerm;
		double gPrime = 2.0 * aTermTheta;
		double h = bTerm * bTerm - 4.0 * aTerm * cTerm;
		double hPrime = 2.0 * bTermTheta * bTerm - 4.0 * aTermTheta * cTerm -
			4.0 * aTerm * cTermTheta;
		double f = -bTerm + sqrt(h);
		double fPrime = -bTermTheta + hPrime / (2.0 * sqrt(h));

		return (fPrime * g - f * gPrime) / (g * g);
	}

	//! Save the total depth
	double totalDepth;

	//! Constants
	const double gasConstant = 8.314; // J K^-1 mol^-1
	const double nAv = 6.022e23;
	const double planckConstant = 6.626e-34; // J s, h
	const double sphereSurface = 12.56e12; // 1 mm radius in nm^2, A_S
	const double purgeFlow = 5.0e20; // nm^3 s^-1, V_p
	const double condensationCoef = 0.4; // sigma
	const double adjSites = 4.0; // N_z
	const double adsorptionE = 15.0e21; // kg nm^2 s^-2 mol^-1
	const double solutionE = 23.0e21; // kg nm^2 s^-2 mol^-1
	const double diffusionE = 141.5e21; // kg nm^2 s^-2 mol^-1
	const double t2MolWeight = 6.032e-3; // kg mol^-1
	const double h2MolWeight = 2.016e-3; // kg mol^-1
	const double htMolWeight = 4.024e-3; // kg mol^-1
	const double surfaceSites = 10.0; // nm^-2, N_S
	const double gamma_H2_H = 2.0;
	const double gamma_T2_T = 2.0;
	const double gamma_HT_H = 1.0;
	const double gamma_HT_T = 1.0;
};

namespace detail
{
class LiReactionGenerator :
	public ReactionGenerator<LiReactionNetwork, LiReactionGenerator>
{
	friend class ReactionGeneratorBase<LiReactionNetwork, LiReactionGenerator>;

public:
	using NetworkType = LiReactionNetwork;
	using Subpaving = typename NetworkType::Subpaving;
	using IndexType = typename NetworkType::IndexType;

	using Superclass =
		ReactionGenerator<LiReactionNetwork, LiReactionGenerator>;

	using Superclass::Superclass;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	operator()(IndexType i, IndexType j, TTag tag) const;

	void
	addConnectivity(Connectivity& conn);

private:
	ReactionCollection<NetworkType>
	getReactionCollection() const;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/LiClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/LiReactionNetwork.tpp>
#endif
