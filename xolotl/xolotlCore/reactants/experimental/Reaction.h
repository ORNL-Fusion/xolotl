#pragma once


namespace xolotlCore
{
namespace experimental
{
template <typename TDerived>
class Reaction
{
public:
    enum class Type
    {
        production,
        dissociation
    };

    void
    computeProductionFluxes(
        Kokkos::View<double*, Kokkos::MemoryUnmanaged> concentrations,
        Kokkos::View<double*, Kokkos::MemoryUnmanaged> fluxes)
    {
        //TODO: compute flux using TDerived definition
        //      auto f = computeProductionFlux(/*_reactants[0], _reacants[1]*/);
        double f = 1.0;
        fluxes[_reactants[0]] = -f;
        fluxes[_reactants[1]] = -f;
        auto invalid = plsm::invalid<std::size_t>;
        for (auto prodId : _products) {
            if (prodId == invalid) {
                continue;
            }
            fluxes[prodId] = f;
        }
    }

    void
    computeDissociationFluxes(
        Kokkos::View<double*, Kokkos::MemoryUnmanaged> concentrations,
        Kokkos::View<double*, Kokkos::MemoryUnmanaged> fluxes)
    {
        //TODO: compute flux using TDerived definition
        //      auto f = computeDissociationFlux(/*_reactants[0]*/);
        double f = 1.0;
        fluxes[_reactants[0]] = -f;
        fluxes[_products[0]] = f;
        fluxes[_products[1]] = f;
    }


    void
    computeFluxes(Kokkos::View<double*, Kokkos::MemoryUnmanaged> concentrations,
        Kokkos::View<double*, Kokkos::MemoryUnmanaged> fluxes)
    {
        //TODO: Use function_ref and make this decision at construction

        switch (_type) {
            case Type::production:
                computeProductionFluxes(concentrations, fluxes);
                break;
            case Type::dissociation:
                computeDissociationFluxes(concentrations, fluxes);
                break;
        }
    }

private:
    Type _type;

    //Cluster indices for LHS and RHS
    //Dissociation reactions always have 1 input and 2 outputs
    //Production reactions always have 2 inputs, but may have 0, 1, or 2 outputs
    //TODO: Both of these should default their elements to some 'invalid' index
    Kokkos::Array<std::size_t, 2> _reactants;
    Kokkos::Array<std::size_t, 2> _products;

    //! Reaction rate, k
    double _rate;
};
}
}
