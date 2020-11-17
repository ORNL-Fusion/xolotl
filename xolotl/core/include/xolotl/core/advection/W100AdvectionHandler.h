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
 * the physical parts for the advection of mobile helium clusters in W100.
 */
class W100AdvectionHandler : public TungstenAdvectionHandler
{
public:
	//! The Constructor
	W100AdvectionHandler() :
		TungstenAdvectionHandler(
			{2.28e-3, 5.06e-3, 7.26e-3, 15.87e-3, 16.95e-3, 27.16e-3, 35.56e-3})
	{
	}

	//! The Destructor
	~W100AdvectionHandler() = default;
};
} /* end namespace advection */
} /* end namespace core */
} /* end namespace xolotl */
