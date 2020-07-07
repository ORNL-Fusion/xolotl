#include <xolotl/factory/solver/SolverFactory.h>

namespace xolotl
{
namespace factory
{
namespace solver
{
SolverFactory::SolverFactory() = default;

SolverFactory::~SolverFactory() = default;

SolverFactory&
SolverFactory::get()
{
	static SolverFactory factory;
	return factory;
}

std::shared_ptr<xolotl::solver::ISolver>
SolverFactory::generateSolver(const options::Options& options)
{
	const auto& solverName = options.getSolverName();
	auto it = _generators.find(solverName);
	if (it == _generators.end()) {
		throw std::runtime_error(
			"No solver found for name \"" + solverName + "\"");
	}
	return it->second(options);
}

bool
SolverFactory::registerGenerator(
	const std::string& name, const Generator& generator)
{
	return _generators.emplace(name, generator).second;
}
} // namespace solver
} // namespace factory
} // namespace xolotl
