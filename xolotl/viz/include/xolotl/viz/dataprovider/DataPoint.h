#pragma once

#include <plsm/SpaceVector.h>

namespace xolotl
{
namespace viz
{
namespace dataprovider
{
struct DataPoint
{
	/**
	 * The time step.
	 */
	double t {};

	/**
	 * The X position on the grid.
	 */
	double x {};

	/**
	 * The Y position on the grid.
	 */
	double y {};

	/**
	 * The Z position on the grid.
	 */
	double z {};

	/**
	 * Scalar field value at time step t and grid position x,y,z.
	 */
	double value {};
};
} // namespace dataprovider
} // namespace viz
} // namespace xolotl
