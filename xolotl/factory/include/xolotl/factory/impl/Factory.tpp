#pragma once

#include <xolotl/factory/Factory.h>
#include <xolotl/options/ConfOptions.h>

namespace xolotl
{
namespace factory
{
template <typename TFactory, typename THandlerBase>
Factory<TFactory, THandlerBase>::Factory()
{
}

template <typename TFactory, typename THandlerBase>
Factory<TFactory, THandlerBase>::~Factory() = default;

template <typename TFactory, typename THandlerBase>
TFactory&
Factory<TFactory, THandlerBase>::get(const UserInitializer& callback)
{
	static TFactory factory;
	callback();
	return factory;
}

template <typename TFactory, typename THandlerBase>
std::shared_ptr<THandlerBase>
Factory<TFactory, THandlerBase>::generate(
	const std::string& name, const options::IOptions& options)
{
	auto it = _generators.find(name);
	if (it == _generators.end()) {
		throw std::runtime_error(TFactory::getFactoryName() +
			": No handler generator found for \"" + name + "\"");
	}
	return it->second(options);
}

template <typename TFactory, typename THandlerBase>
std::shared_ptr<THandlerBase>
Factory<TFactory, THandlerBase>::generate(const std::string& name)
{
	return generate(name, options::ConfOptions{});
}

template <typename TFactory, typename THandlerBase>
std::shared_ptr<THandlerBase>
Factory<TFactory, THandlerBase>::generate(const options::IOptions& options)
{
	return generate(TFactory::getName(options), options);
}

template <typename TFactory, typename THandlerBase>
bool
Factory<TFactory, THandlerBase>::registerGenerator(
	const std::string& name, const Generator& generator)
{
	return _generators.emplace(name, generator).second;
}
} // namespace factory
} // namespace xolotl
