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
namespace temperature
{
class ITemperatureHandler;
}
} // namespace core

namespace factory
{
namespace temperature
{
class TemperatureHandlerFactory
{
public:
	using Generator =
		std::function<std::shared_ptr<core::temperature::ITemperatureHandler>(
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

	TemperatureHandlerFactory(const TemperatureHandlerFactory&) = delete;

	static TemperatureHandlerFactory&
	get();

	std::shared_ptr<core::temperature::ITemperatureHandler>
	generateTemperatureHandler(const options::Options& options);

	bool
	registerGenerator(const std::string& name, const Generator& generator);

private:
	TemperatureHandlerFactory();

	~TemperatureHandlerFactory();

private:
	std::unordered_map<std::string, Generator> _generators;
};

template <typename THandler>
TemperatureHandlerFactory::Registration<THandler>::Registration(
	const std::string& name)
{
	TemperatureHandlerFactory::get().registerGenerator(
		name, [](const options::Options& options) {
			return std::make_shared<THandler>(options);
		});
}

template <typename THandler>
TemperatureHandlerFactory::RegistrationCollection<
	THandler>::RegistrationCollection(const std::vector<std::string>& names)
{
	for (const auto& name : names) {
		registrations.emplace_back(name);
	}
}
} // end namespace temperature
} // end namespace factory
} // end namespace xolotl
