#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include <xolotl/options/Options.h>

namespace xolotl
{
namespace solver
{
class ISolver;
}

namespace factory
{
namespace solver
{
class SolverFactory
{
public:
	using Generator = std::function<std::shared_ptr<xolotl::solver::ISolver>(
		const options::Options& options)>;

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

	SolverFactory(const SolverFactory&) = delete;

	static SolverFactory&
	get();

	std::shared_ptr<xolotl::solver::ISolver>
	generateSolver(const options::Options& options);

	bool
	registerGenerator(const std::string& name, const Generator& generator);

private:
	SolverFactory();

	~SolverFactory();

private:
	std::unordered_map<std::string, Generator> _generators;
};

template <typename THandler>
SolverFactory::Registration<THandler>::Registration(const std::string& name)
{
	SolverFactory::get().registerGenerator(
		name, [](const options::Options& options) {
			return std::make_shared<THandler>(options);
		});
}

template <typename THandler>
SolverFactory::RegistrationCollection<THandler>::RegistrationCollection(
	const std::vector<std::string>& names)
{
	for (const auto& name : names) {
		registrations.emplace_back(name);
	}
}
} // namespace solver
} // namespace factory
} // namespace xolotl
