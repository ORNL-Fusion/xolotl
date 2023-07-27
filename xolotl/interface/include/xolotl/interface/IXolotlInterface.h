#pragma once

namespace xolotl
{
namespace interface
{
class IXolotlInterface
{
public:
	virtual ~IXolotlInterface()
	{
	}

	/**
	 * Run the PETSc solve
	 */
	virtual void
	solveXolotl() = 0;
};
} // namespace interface
} // namespace xolotl
