#include <xolotl/interface/ComputeContext.h>
#include <xolotl/interface/Interface.h>
#include <xolotl/interface/MultiXolotl.h>
#include <xolotl/interface/XolotlInterface.h>
#include <xolotl/options/Options.h>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace interface
{
std::unique_ptr<IXolotlInterface>
makeXolotlInterface(int& argc, const char* argv[])
{
	// Init HPC context
	auto context = std::make_shared<ComputeContext>(argc, argv);
	util::setMPIComm(MPI_COMM_WORLD);

	// Get options
	auto options = std::make_shared<options::Options>();
	options->readParams(argc, argv);

	// Create instance(s)
	if (options->useSubnetworks()) {
		return std::make_unique<MultiXolotl>(context, options);
	}
	else {
		return std::make_unique<XolotlInterface>(context, options);
	}
}
} // namespace interface
} // namespace xolotl
