#include <xolotl/core/network/PSINetworkHandler.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template auto
generatePSIReactionNetwork<PSIFullSpeciesList>(const options::Options& options);
}

template class PSINetworkHandler<PSIFullSpeciesList>;

namespace detail
{
auto psiFullSpeciesListNetworkHandlerRegistrations =
	xolotl::factory::network::NetworkHandlerFactory::RegistrationCollection<
		PSINetworkHandler<PSIFullSpeciesList>>(
		{"W100", "W110", "W111", "W211", "TRIDYN", "Pulsed"});
}
} // namespace network
} // namespace core
} // namespace xolotl
