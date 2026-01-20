#include <xolotl/factory/interface/MaterialSubOptionsFactory.h>
#include <xolotl/interface/AlphaZrSubOptions.h>

namespace xolotl
{
namespace interface
{
auto alphaZrSubOptionsRegistrations =
	xolotl::factory::interface::MaterialSubOptionsFactory::
		RegistrationCollection<AlphaZrSubOptions>({"AlphaZr"});

AlphaZrSubOptions::AlphaZrSubOptions(const options::IOptions& options) :
	_options(options)
{
}

std::vector<std::shared_ptr<options::IOptions>>
AlphaZrSubOptions::getSubOptions() const
{
	std::vector<std::shared_ptr<options::IOptions>> ret;
	const auto& netParams = _options.getNetworkParameters();
	auto tmpParams = netParams;
	for (std::size_t i = 0; i < netParams.size(); ++i) {
		if (netParams[i] == 0) {
			continue;
		}

		tmpParams.assign(netParams.size(), 0);
		tmpParams[i] = netParams[i];
		auto subOpts = _options.makeCopy();
		subOpts->setNetworkParameters(tmpParams);

		ret.push_back(subOpts);
	}
	return ret;
}
} // namespace interface
} // namespace xolotl
