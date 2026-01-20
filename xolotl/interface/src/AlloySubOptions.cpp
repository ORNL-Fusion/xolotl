#include <xolotl/factory/interface/MaterialSubOptionsFactory.h>
#include <xolotl/interface/AlloySubOptions.h>

namespace xolotl
{
namespace interface
{
auto alloySubOptionsRegistrations =
	xolotl::factory::interface::MaterialSubOptionsFactory::
		RegistrationCollection<AlloySubOptions>({"800H"});

AlloySubOptions::AlloySubOptions(const options::IOptions& options) :
	_options(options)
{
}

std::vector<std::shared_ptr<options::IOptions>>
AlloySubOptions::getSubOptions() const
{
	std::vector<std::shared_ptr<options::IOptions>> ret;
	const auto& netParams = _options.getNetworkParameters();

	auto tmpParams = netParams;

	// Vacancy
	if (netParams[3] != 0) {
		tmpParams[0] = netParams[0];
		tmpParams[1] = netParams[1];
		tmpParams[2] = 0;
		tmpParams[3] = netParams[3];
		tmpParams[4] = 0;
		ret.emplace_back(_options.makeCopy())->setNetworkParameters(tmpParams);
	}

	// Interstitial
	if (netParams[4] != 0) {
		tmpParams[0] = netParams[0];
		tmpParams[1] = 0;
		tmpParams[2] = 0;
		tmpParams[3] = 0;
		tmpParams[4] = netParams[4];
		ret.emplace_back(_options.makeCopy())->setNetworkParameters(tmpParams);
	}

	return ret;
}
} // namespace interface
} // namespace xolotl
