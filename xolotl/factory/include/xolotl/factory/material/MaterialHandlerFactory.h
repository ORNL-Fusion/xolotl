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
namespace material
{
class IMaterialHandler;
}
} // namespace core

namespace factory
{
namespace material
{
class MaterialHandlerFactory
{
public:
	using Generator =
		std::function<std::shared_ptr<core::material::IMaterialHandler>(
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

	MaterialHandlerFactory(const MaterialHandlerFactory&) = delete;

	static MaterialHandlerFactory&
	get();

	std::shared_ptr<core::material::IMaterialHandler>
	generateMaterialHandler(const options::Options& options);

	bool
	registerGenerator(const std::string& name, const Generator& generator);

private:
	MaterialHandlerFactory();

	~MaterialHandlerFactory();

private:
	std::unordered_map<std::string, Generator> _generators;
};

template <typename THandler>
MaterialHandlerFactory::Registration<THandler>::Registration(
	const std::string& name)
{
	MaterialHandlerFactory::get().registerGenerator(
		name, [](const options::Options& options) {
			return std::make_shared<THandler>(options);
		});
}

template <typename THandler>
MaterialHandlerFactory::RegistrationCollection<
	THandler>::RegistrationCollection(const std::vector<std::string>& names)
{
	for (const auto& name : names) {
		registrations.emplace_back(name);
	}
}
} // namespace material
} // namespace factory
} // namespace xolotl
