#ifndef VIZHANDLERREGISTRYFACTORY_H
#define VIZHANDLERREGISTRYFACTORY_H

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <xolotl/options/Options.h>

namespace xolotl
{
namespace viz
{
class IVizHandlerRegistry;
}

namespace factory
{
namespace viz
{
class VizHandlerRegistryFactory
{
public:
	using Generator =
		std::function<std::shared_ptr<xolotl::viz::IVizHandlerRegistry>(
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

	VizHandlerRegistryFactory(const VizHandlerRegistryFactory&) = delete;

	static VizHandlerRegistryFactory&
	get();

	std::shared_ptr<xolotl::viz::IVizHandlerRegistry>
	generateVizHandlerRegistry(const options::Options& options);

	bool
	registerGenerator(const std::string& name, const Generator& generator);

private:
	VizHandlerRegistryFactory();

	~VizHandlerRegistryFactory();

private:
	std::unordered_map<std::string, Generator> _generators;
};

template <typename THandler>
VizHandlerRegistryFactory::Registration<THandler>::Registration(
	const std::string& name)
{
	VizHandlerRegistryFactory::get().registerGenerator(
		name, [](const options::Options& options) {
			return std::make_shared<THandler>(options);
		});
}

template <typename THandler>
VizHandlerRegistryFactory::RegistrationCollection<
	THandler>::RegistrationCollection(const std::vector<std::string>& names)
{
	for (const auto& name : names) {
		registrations.emplace_back(name);
	}
}
} // end namespace viz
} // end namespace factory
} // end namespace xolotl

#endif // VIZHANDLERREGISTRYFACTORY_H
