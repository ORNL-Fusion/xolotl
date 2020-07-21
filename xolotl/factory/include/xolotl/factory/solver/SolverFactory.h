#pragma once

#include <xolotl/factory/Factory.h>

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
class SolverFactory : public Factory<SolverFactory, xolotl::solver::ISolver>
{
public:
	static std::string
	getFactoryName() noexcept
	{
		return "SolverFactory";
	}

	static std::string
	getName(const options::Options& options)
	{
		return options.getSolverName();
	}
};
} // namespace solver
} // namespace factory
} // namespace xolotl
