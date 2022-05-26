#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace factory
{
template <typename TFactory, typename THandlerBase>
class Factory
{
public:
	using Generator =
		std::function<std::shared_ptr<THandlerBase>(const options::IOptions&)>;

	using UserInitializer = std::function<void(void)>;

	template <typename THandler>
	static auto
	makeDefaultGenerator() noexcept
	{
		return [](const options::IOptions& options) {
			return std::make_shared<THandler>(options);
		};
	}

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

	Factory(const Factory&) = delete;

	static TFactory&
	get(const UserInitializer& callback = []() {});

	std::shared_ptr<THandlerBase>
	generate(const std::string& name, const options::IOptions& options);

	std::shared_ptr<THandlerBase>
	generate(const std::string& name);

	std::shared_ptr<THandlerBase>
	generate(const options::IOptions& options);

	bool
	registerGenerator(const std::string& name, const Generator& generator);

protected:
	Factory();

	~Factory();

private:
	std::unordered_map<std::string, Generator> _generators;
};

template <typename TFactory, typename THandlerBase>
template <typename THandler>
Factory<TFactory, THandlerBase>::Registration<THandler>::Registration(
	const std::string& name)
{
	TFactory::get().registerGenerator(
		name, TFactory::template makeDefaultGenerator<THandler>());
}

template <typename TFactory, typename THandlerBase>
template <typename THandler>
Factory<TFactory, THandlerBase>::RegistrationCollection<
	THandler>::RegistrationCollection(const std::vector<std::string>& names)
{
	for (const auto& name : names) {
		registrations.emplace_back(name);
	}
}
} // namespace factory
} // namespace xolotl
