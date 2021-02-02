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
protected:
	//! Collection of diffusing clusters.
	std::vector<std::size_t> diffusingClusters;

	/**
	 * The number of degrees of freedom in the network
	 */
	int dof;

	/**
	 * The x grid
	 */
	std::vector<double> xGrid;

	/**
	 * The heat flux in W.m-2
	 */
	double heatFlux;

	/**
	 * The heat conductivity
	 */
	double heatConductivity;

	/**
	 * The surface position
	 */
	int surfacePosition;

	/**
	 * The beta factor
	 */
	double beta;

	/**
	 * The He flux (nm-2 s-1)
	 */
	double J;

public:
	//! The Constructor
	DummySoretDiffusionHandler() :
		dof(0),
		surfacePosition(0),
		heatConductivity(0.0),
		heatFlux(0.0),
		beta(0.0065),
		J(1.0e3)
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
