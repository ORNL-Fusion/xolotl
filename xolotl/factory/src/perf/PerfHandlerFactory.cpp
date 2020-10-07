#include <xolotl/factory/impl/Factory.tpp>
#include <xolotl/factory/perf/PerfHandlerFactory.h>

namespace xolotl
{
namespace factory
{
using perf::PerfHandlerFactory;
using xolotl::perf::IPerfHandler;

template PerfHandlerFactory&
Factory<PerfHandlerFactory, IPerfHandler>::get();

template std::shared_ptr<IPerfHandler>
Factory<PerfHandlerFactory, IPerfHandler>::generate(const options::Options&);

template bool
Factory<PerfHandlerFactory, IPerfHandler>::registerGenerator(
	const std::string&, const Generator&);
} // namespace factory
} // namespace xolotl
