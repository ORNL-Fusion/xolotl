#pragma once

#include <plsm/EnumIndexed.h>

#include <xolotl/core/Constants.h>
#include <xolotl/core/MathUtils.h>

namespace xolotl
{
namespace core
{
namespace network
{
template <typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() >= 2),
    double >
getRateForSink(const double r, const double dc, const double latticeParameter)
{
    double strength = ::xolotl::core::alloysinkStrength * dc;

    return strength;
}

template <typename TRegion>
KOKKOS_INLINE_FUNCTION
std::enable_if_t<(numberOfVacancySpecies<typename TRegion::EnumIndex>() < 2),
    double >
getRateForSink(const double r, const double dc, const double latticeParameter)
{
    double r0 = latticeParameter * 0.75 * sqrt(3.0);
    double rho = 0.0003;
    constexpr double pi = ::xolotl::core::pi;

    double strength = -4.0 * pi * rho
            / log(pi * rho * (r + r0)
            * (r + r0)) * dc;

    return strength;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
double
SinkReaction<TNetwork, TDerived>::computeRate(IndexType gridIndex)
{
    auto cl = this->_clusterData.getCluster(_reactant);
    double r = cl.getReactionRadius();
    double dc = cl.getDiffusionCoefficient(gridIndex);
    using Types = detail::ReactionNetworkTypes<TNetwork>;
    using Region = typename Types::Region;

    double strength = this->asDerived()->getSinkBias() * getRateForSink<Region>(r, dc, this->_clusterData.getLatticeParameter());
    
    return strength;
}
}
}
}
