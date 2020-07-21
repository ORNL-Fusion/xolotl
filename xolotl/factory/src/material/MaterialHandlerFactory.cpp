#include <xolotl/factory/impl/Factory.tpp>
#include <xolotl/factory/material/MaterialHandlerFactory.h>

namespace xolotl
{
namespace factory
{
using core::material::IMaterialHandler;
using material::MaterialHandlerFactory;

template MaterialHandlerFactory&
Factory<MaterialHandlerFactory, IMaterialHandler>::get();

template std::shared_ptr<IMaterialHandler>
Factory<MaterialHandlerFactory, IMaterialHandler>::generate(
	const options::Options&);

template bool
Factory<MaterialHandlerFactory, IMaterialHandler>::registerGenerator(
	const std::string&, const Generator&);
} // namespace factory
} // namespace xolotl
