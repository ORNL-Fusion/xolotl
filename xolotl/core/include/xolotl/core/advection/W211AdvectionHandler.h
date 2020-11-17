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
 * the physical parts for the advection of mobile helium clusters in W211.
 */
class W211AdvectionHandler : public TungstenAdvectionHandler
{
public:
	//! The Constructor
	W211AdvectionHandler() :
		TungstenAdvectionHandler({1.49e-3, 3.69e-3, 12.34e-3, 14.11e-3,
			19.14e-3, 35.77e-3, 67.65e-3})
	{
	}

	//! The Destructor
	~W211AdvectionHandler() = default;
};
} /* end namespace advection */
} /* end namespace core */
} /* end namespace xolotl */
