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
	virtual void
	initialize(const network::IReactionNetwork& network,
		network::IReactionNetwork::SparseFillMap& ofill,
		network::IReactionNetwork::SparseFillMap& dfill,
		std::vector<double> grid) override
	{
		return;
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	void
	setHeatFlux(double flux) override
	{
		return;
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	void
	setHeatConductivity(double cond) override
	{
		return;
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
	setInterfaceLocation(double loc) override
	{
		return;
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	virtual void
	computeDiffusion(network::IReactionNetwork& network, double** concVector,
		double* updatedConcOffset, double hxLeft, double hxRight, int ix,
		double sy = 0.0, int iy = 0, double sz = 0.0, int iz = 0) const
	{
		return;
	}

	/**
	 * \see ISoretDiffusionHandler.h
	 */
	virtual bool
	computePartialsForDiffusion(network::IReactionNetwork& network,
		double** concVector, double* val, int* indices, double hxLeft,
		double hxRight, int ix, double sy = 0.0, int iy = 0, double sz = 0.0,
		int iz = 0) const
	{
		return false;
	}
};
// end class DummySoretDiffusionHandler

} /* end namespace modified */
} /* end namespace core */
} /* end namespace xolotl */
#endif
