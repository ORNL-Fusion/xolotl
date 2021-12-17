#include <xolotl/core/temperature/ConstantHandler.h>
#include <xolotl/factory/temperature/TemperatureHandlerFactory.h>
#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
namespace detail
{
auto constantTemperatureHandlerRegistration =
	xolotl::factory::temperature::TemperatureHandlerFactory::Registration<
		ConstantHandler>("constant");
}

ConstantHandler::ConstantHandler(double constTemperature) :
	temperature(constTemperature)
{
	auto xolotlComm = util::getMPIComm();
	int procId;
	MPI_Comm_rank(xolotlComm, &procId);
	if (procId == 0) {
		XOLOTL_LOG << "TemperatureHandler: Using a constant temperature of: "
				   << temperature << " K";
	}
}

ConstantHandler::ConstantHandler(const options::IOptions& options) :
	ConstantHandler(options.getTempParam())
{
}

ConstantHandler::~ConstantHandler()
{
}
} // namespace temperature
} // namespace core
} // namespace xolotl
