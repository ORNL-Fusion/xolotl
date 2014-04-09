#ifndef TCLUSTER_H
#define TCLUSTER_H

// Includes
#include "HCluster.h"

namespace xolotlCore {

//! This class represents a cluster composed entirely of tritium.
class TCluster : public HCluster {

private:

public:

	//! The constructor. All TClusters must be initialized with a size.
	TCluster(int nT, std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	//! Destructor
	~TCluster();

};
//end class TCluster

} /* end namespace xolotlCore */
#endif
