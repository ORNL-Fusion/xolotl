#pragma once

#include <xolotl/core/network/FeReaction.h>
#include <xolotl/core/network/FeTraits.h>
#include <xolotl/core/network/ReactionNetwork.h>
#include <xolotl/util/MathUtils.h>
#include <xolotl/core/Constants.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
class FeReactionGenerator;
}

class FeReactionNetwork : public ReactionNetwork<FeReactionNetwork>
{
	friend class ReactionNetwork<FeReactionNetwork>;

public:
	using Superclass = ReactionNetwork<FeReactionNetwork>;
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

	void setReactionParams(std::string) override;

	std::string
	getMonitorOutputFileName() const override
	{
		return "Fe.dat";
	}

	void
        initializeExtraClusterData(const options::IOptions& options);

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
        computeBubbleRadius(
                double amount, double latticeParameter, IndexType typeSwitch)
        {
                // Find the edge of the phase space
                const auto& largestReg = this->getCluster(largestClusterId).getRegion();
                Composition hiLargest = largestReg.getUpperLimitPoint();
                double largestSize =100.0;// hiLargest[Species::V] - 1;

                // Get the minimum amount for a valid value
                amount = util::max(amount, largestSize);
        	double radius = 0.0;

		return latticeParameter * pow((3.0 * amount) / ::xolotl::core::pi, (1.0 / 3.0)) * 0.5; // matches V in Reaction bubble radius

		// Idea for section to match bubble radius with how reaction radius is found in reaction bubble radius, current implimentation
		// yields unphysical results, but could edit this later. For now it is all skipped and defaults to V radius calcuation. More so 
		// for Std model than SSBM
	 	switch (typeSwitch) {
		// Trap
		case 0:
			{
				radius = 0.0;
			}
		// Heliem	
		case 1:
			{
				double FourPi = 4.0 * ::xolotl::core::pi;
        	                double aCubed = pow(latticeParameter, 3);
	                        double termOne =
                                	pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed * amount,
                        	                (1.0 / 3.0));
                	        double termTwo =
        	                        pow((3.0 / FourPi) * (1.0 / 10.0) * aCubed, (1.0 / 3.0));
	                        radius = heliumRadius + termOne - termTwo;
			}
		// Bubble (V)
		case 2:
			{
				radius = latticeParameter *
                                	pow((3.0 * amount) / ::xolotl::core::pi,
                                        	(1.0 / 3.0)) * 0.5;
			}
		// Interstitial
		case 3:
			{
				radius =  latticeParameter * pow(3.0 / ::xolotl::core::pi, (1.0 / 3.0)) * 0.5;
			}
		}       
        	
		return radius;
	}

	std::string
	getMonitorDataHeaderString() const override;

	void
	addMonitorDataValues(Kokkos::View<const double*> conc, double fac,
		std::vector<double>& totalVals) override;

	std::size_t
	getMonitorDataLineSize() const override
	{
		return 2 + getSpeciesListSize() * 12;
	}

	void
	writeMonitorDataLine(
		const std::vector<double>& localData, double time) override;

public:
	IndexType largestClusterId;

private:
	double
	checkLatticeParameter(double latticeParameter);

	double
	computeAtomicVolume(double latticeParameter)
	{
		// 2 atoms per cell
		return 0.5 * latticeParameter * latticeParameter * latticeParameter;
	}

	double
	checkImpurityRadius(double impurityRadius);

	detail::FeReactionGenerator
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
};

namespace detail
{
class FeReactionGenerator :
	public ReactionGenerator<FeReactionNetwork, FeReactionGenerator>
{
	friend class ReactionGeneratorBase<FeReactionNetwork, FeReactionGenerator>;

public:
	using NetworkType = FeReactionNetwork;
	using Subpaving = typename NetworkType::Subpaving;
	using IndexType = typename NetworkType::IndexType;

	using Superclass =
		ReactionGenerator<FeReactionNetwork, FeReactionGenerator>;

	using Superclass::Superclass;

	FeReactionGenerator(const FeReactionNetwork& network) :
		Superclass(network),
		largestClusterId(network.largestClusterId)
	{
	}

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	operator()(IndexType i, IndexType j, TTag tag) const;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addSinks(IndexType i, TTag tag) const;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addSingleSizeReactions(IndexType i, IndexType j, TTag tag) const;

	template <typename TTag>
	KOKKOS_INLINE_FUNCTION
	void
	addTraps(IndexType i, IndexType j, TTag tag) const;

private:
	ReactionCollection<NetworkType>
	getReactionCollection() const;
	
	IndexType largestClusterId;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl

#include <xolotl/core/network/FeClusterGenerator.h>

#if defined(XOLOTL_INCLUDE_RN_TPP_FILES)
#include <xolotl/core/network/impl/FeReactionNetwork.tpp>
#endif
