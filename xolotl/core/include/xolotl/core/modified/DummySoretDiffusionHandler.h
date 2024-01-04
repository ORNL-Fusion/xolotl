#ifndef DUMMYSORETDIFFUSIONHANDLER_H
#define DUMMYSORETDIFFUSIONHANDLER_H

// Includes
#include <xolotl/core/modified/ISoretDiffusionHandler.h>

namespace xolotl
{
namespace core
{
namespace modified
{
/**
 * This class realizes the IDiffusionHandler interface responsible for all
 * the physical parts for the diffusion of mobile clusters. This class doesn't
 * do anything.
 */
class DummySoretDiffusionHandler : public ISoretDiffusionHandler
{
public:
	//! The Constructor
	DummySoretDiffusionHandler()
	{
	}

	//! The Destructor
	~DummySoretDiffusionHandler()
	{
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	void
	initialize(network::IReactionNetwork& network,
		std::vector<core::RowColPair>& idPairs, std::vector<double> grid,
		int xs) override
	{
		return;
	}

	int
	getNumberOfDiffusing() const override
	{
		return 0;
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	void
	updateSurfacePosition(int surfacePos) override
	{
		return;
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	void
	computeDiffusion(network::IReactionNetwork& network,
		const StencilConcArray& concVector,
		Kokkos::View<double*> updatedConcOffset, double hxLeft, double hxRight,
		int ix, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) const override
	{
		return;
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	bool
	computePartialsForDiffusion(network::IReactionNetwork& network,
		const StencilConcArray& concVector, Kokkos::View<double*> values,
		double hxLeft, double hxRight, int ix, double sy = 0.0, int iy = 0,
		double sz = 0.0, int iz = 0) const override
	{
		return false;
	}
};
// end class DummySoretDiffusionHandler

} /* end namespace modified */
} /* end namespace core */
} /* end namespace xolotl */
#endif
