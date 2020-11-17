#pragma once

#include <xolotl/core/advection/TungstenAdvectionHandler.h>

namespace xolotl
{
namespace core
{
namespace advection
{
/**
 * This class realizes the IAdvectionHandler interface responsible for all
 * the physical parts for the advection of mobile helium cluster.
 */
class W110AdvectionHandler : public TungstenAdvectionHandler
{
public:
	//! The Constructor
	W110AdvectionHandler() :
		TungstenAdvectionHandler(
			{0.92e-3, 1.48e-3, 6.73e-3, 6.18e-3, 33.61e-3, 37.58e-3, 41.90e-3})
	{
	}

	//! The Destructor
	~W110AdvectionHandler() = default;
};
} /* end namespace advection */
} /* end namespace core */
} /* end namespace xolotl */
