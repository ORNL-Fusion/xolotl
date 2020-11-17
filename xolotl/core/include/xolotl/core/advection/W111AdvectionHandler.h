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
 * the physical parts for the advection of mobile helium clusters in W111.
 */
class W111AdvectionHandler : public TungstenAdvectionHandler
{
public:
	//! The Constructor
	W111AdvectionHandler() :
		TungstenAdvectionHandler(
			{3.65e-3, 6.40e-3, 16.38e-3, 9.84e-3, 44.40e-3, 52.12e-3, 81.57e-3})
	{
	}

	//! The Destructor
	~W111AdvectionHandler() = default;
};
} /* end namespace advection */
} /* end namespace core */
} /* end namespace xolotl */
