#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <xolotl/options/Options.h>

namespace xolotl
{
namespace core
{
namespace network
{
class INetworkHandler;
}
} // namespace core

namespace factory
{
namespace network
{
class NetworkHandlerFactory
{
public:
	using Generator =
		std::function<std::shared_ptr<core::network::INetworkHandler>(
			const options::Options&)>;

	template <typename THandler>
	struct Registration
	{
		Registration(const std::string& name);
	};

	template <typename THandler>
	struct RegistrationCollection
	{
		RegistrationCollection(const std::vector<std::string>& names);

		std::vector<Registration<THandler>> registrations;
	};

	NetworkHandlerFactory(const NetworkHandlerFactory&) = delete;

	static NetworkHandlerFactory&
	get();

	std::shared_ptr<core::network::INetworkHandler>
	generateNetworkHandler(const options::Options& options);

	bool
	registerGenerator(const std::string& name, const Generator& generator);

private:
	NetworkHandlerFactory();

	~NetworkHandlerFactory();

private:
	std::unordered_map<std::string, Generator> _generators;
};

template <typename THandler>
NetworkHandlerFactory::Registration<THandler>::Registration(
	const std::string& name)
{
	NetworkHandlerFactory::get().registerGenerator(
		name, [](const options::Options& options) {
			return std::make_shared<THandler>(options);
		});
}

template <typename THandler>
NetworkHandlerFactory::RegistrationCollection<THandler>::RegistrationCollection(
	const std::vector<std::string>& names)
{
	for (const auto& name : names) {
		registrations.emplace_back(name);
	}
}
} // namespace network
} // namespace factory
} // namespace xolotl
