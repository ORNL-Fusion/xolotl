#include <xolotl/factory/impl/Factory.tpp>
#include <xolotl/factory/solver/SolverFactory.h>

namespace xolotl
{
namespace factory
{
using solver::SolverFactory;
using xolotl::solver::ISolver;

template SolverFactory&
Factory<SolverFactory, ISolver>::get();

template std::shared_ptr<ISolver>
Factory<SolverFactory, ISolver>::generate(const options::Options&);

template bool
Factory<SolverFactory, ISolver>::registerGenerator(
	const std::string&, const Generator&);
} // namespace factory
} // namespace xolotl
