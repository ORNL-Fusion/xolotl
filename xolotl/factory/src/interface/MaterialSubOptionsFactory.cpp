#include <xolotl/factory/impl/Factory.tpp>
#include <xolotl/factory/interface/MaterialSubOptionsFactory.h>

namespace xolotl
{
namespace factory
{
using ::xolotl::interface::IMaterialSubOptions;
using ::xolotl::factory::interface::MaterialSubOptionsFactory;

template MaterialSubOptionsFactory&
Factory<MaterialSubOptionsFactory, IMaterialSubOptions>::get(
	const UserInitializer&);

template std::shared_ptr<IMaterialSubOptions>
Factory<MaterialSubOptionsFactory, IMaterialSubOptions>::generate(
	const options::IOptions&);

template bool
Factory<MaterialSubOptionsFactory, IMaterialSubOptions>::registerGenerator(
	const std::string&, const Generator&);
} // namespace factory
} // namespace xolotl
