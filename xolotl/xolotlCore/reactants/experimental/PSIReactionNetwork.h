#pragma once

#include <experimental/ReactionNetwork.h>

namespace xolotlCore
{
namespace experimental
{
template <typename TSpeciesEnum>
class PSIReactionNetwork;


template <typename TSpeciesEnum>
class PSIReaction;


enum class PSIFullSpeciesList
{
    He,
    D,
    T,
    V,
    I
};


template <>
struct hasInterstitial<PSIFullSpeciesList> : std::true_type { };


template <typename TSpeciesEnum>
struct ReactionNetworkTraits<PSIReactionNetwork<TSpeciesEnum>>
{
    using Species = TSpeciesEnum;

    static constexpr std::size_t numSpecies = 5;

    using ReactionType = PSIReaction<TSpeciesEnum>;
};


template <typename TSpeciesEnum>
class PSIClusterGenerator;


template <typename TSpeciesEnum>
class PSIReactionNetwork :
    public ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>
{
public:
    using Superclass = ReactionNetwork<PSIReactionNetwork<TSpeciesEnum>>;
    using Subpaving = typename Superclass::Subpaving;

    using Superclass::Superclass;

    PSIReactionNetwork(Subpaving&& subpaving, std::size_t gridSize,
            const IOptions& options)
        :
        Superclass(std::move(subpaving), gridSize, options,
            PSIClusterGenerator<TSpeciesEnum>{options})
    {
        if (this->getLatticeParameter() <= 0.0) {
            this->setLatticeParameter(tungstenLatticeConstant);
        }

        if (this->getImpurityRadius() <= 0.0) {
            this->setImpurityRadius(heliumRadius);
        }
    }
};


template <typename TSpeciesEnum>
class PSIClusterGenerator
{
public:
    using Species = TSpeciesEnum;
    using NetworkType = PSIReactionNetwork<Species>;
    using Cluster = typename NetworkType::Cluster;
    using Composition = typename NetworkType::Composition;

    PSIClusterGenerator(const IOptions& options)
        :
        _hydrogenRadiusFactor(options.getHydrogenFactor())
    {
    }

    double
    getFormationEnergy(const Cluster& cluster) const noexcept
    {
        const auto& reg = cluster.getRegion();
        if (reg.isSimplex()) {
            Composition comp = reg.getOrigin();
            if (comp.isOnAxis(Species::I)) {
                if (comp[Species::I] < _iFormationEnergies.size()) {
                }
            }
        }
        return 0.0;
    }

    double
    getMigrationEnergy(const Cluster& cluster) const noexcept
    {
        return 0.0;
    }

    double
    getDiffusionFactor(const Cluster& cluster) const noexcept
    {
        return 0.0;
    }

private:
    // I formation energies in eV
    static constexpr Kokkos::Array<double, 6> _iFormationEnergies = {
        10.0, 18.5, 27.0, 35.0, 42.5, 48.0
    };
    // I diffusion factors in nm^2/s
    static constexpr double iDiffusion[] = {
        8.8e+10, 8.0e+10, 3.9e+10, 2.0e+10, 1.0e+10
    };
    // I migration energies in eV
    static constexpr double iMigration[] = { 0.01, 0.02, 0.03, 0.04, 0.05 };

    // He formation energies in eV
    static constexpr double heFormationEnergies[] = {
        6.15, 11.44, 16.35, 21.0, 26.1, 30.24, 34.93, 38.80
    };
    // He diffusion factors in nm^2/s
    static constexpr double heDiffusion[] = {
        2.9e+10, 3.2e+10, 2.3e+10, 1.7e+10, 5.0e+09, 1.0e+09, 5.0e+08
    };
    // He migration energies in eV
    static constexpr double heMigration[] = {
        0.13, 0.20, 0.25, 0.20, 0.12, 0.3, 0.4
    };

    // The diffusion factor for a single deuterium.
    static constexpr double dOneDiffusionFactor = 2.83e+11;
    // The migration energy for a single deuterium.
    static constexpr double dOneMigrationEnergy = 0.38;

    // The diffusion factor for a single tritium.
    static constexpr double tOneDiffusionFactor = 2.31e+11;
    // The migration energy for a single tritium.
    static constexpr double tOneMigrationEnergy = 0.38;

    // The diffusion factor for a single vacancy in nm^2/s
    static constexpr double vOneDiffusion = 1.8e+12;
    // The migration energy for a single vacancy in eV
    static constexpr double vOneMigration = 1.30;

    // The factor between He and H radius sizes
    double _hydrogenRadiusFactor = 0.25;

    /**
     * The maximum number of helium atoms that can be combined with a
     * vacancy cluster with size equal to the index i in the array plus one.
     * For example, an HeV size cluster with size 1 would have
     * size = i+1 = 1 and i = 0. It could support a mixture of up to nine
     * helium atoms with one vacancy.
     */
    static constexpr int maxHePerV[] = {
        9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85,
        90, 95, 98, 100, 101, 103, 105, 107, 109, 110, 112, 116
    };
};
}
}

#include <experimental/PSIReaction.h>
