#include <xolotl/factory/impl/Factory.tpp>
#include <xolotl/factory/perf/PerfHandlerFactory.h>

namespace xolotl
{
namespace factory
{
using perf::PerfHandlerFactory;
using xolotl::perf::IPerfHandler;

template PerfHandlerFactory&
Factory<PerfHandlerFactory, IPerfHandler>::get(const UserInitializer&);

template std::shared_ptr<IPerfHandler>
Factory<PerfHandlerFactory, IPerfHandler>::generate(const options::IOptions&);

template std::shared_ptr<IPerfHandler>
Factory<PerfHandlerFactory, IPerfHandler>::generate(const std::string&);

template bool
Factory<PerfHandlerFactory, IPerfHandler>::registerGenerator(
	const std::string&, const Generator&);
} // namespace factory
} // namespace xolotl
