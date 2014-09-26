#ifndef SIMPLEREACTIONNETWORK_H_
#define SIMPLEREACTIONNETWORK_H_

#include <PSIClusterReactionNetwork.h>
#include <xolotlPerf.h>
#include <DummyHandlerRegistry.h>

namespace testUtils {

/**
 * This class creates a simple reaction network used for testing. It contains
 * 10 each of He, V and I clusters and twenty HeV clusters. The HeV clusters are
 * stored in ascending order of He and then V, (1,1;2,1;2,2;3,2;3,3;etc.).
 *
 * It does not register itself as the ReactionNetwork for its clusters because
 * of limitations with shared_ptrs and "this." So, the
 * TestUtils::getSimpleReactionNetwork() operation should always be called to
 * insure that it is properly initialized.
 */
class SimpleReactionNetwork: public xolotlCore::PSIClusterReactionNetwork {

public:
	/**
	 * Constructor
	 * @param maxClusterSize the maximal size of the clusters that will be in
	 * the network. Set to 10 by default.
	 */
	SimpleReactionNetwork(const int maxClusterSize = 10,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry = std::make_shared<xolotlPerf::DummyHandlerRegistry>());

	//! Destructor
	virtual ~SimpleReactionNetwork();

};

/**
 * This operation creates a SimpleReactionNetwork and makes sure that it is
 * properly registered with the clusters it contains. This operation should
 * always be called instead of constructing a SimpleReactionNetwork manually.
 * @param maxClusterSize the maximal size of the clusters that will be in
 * the network. Set to 10 by default.
 * @return The reaction network.
 */
std::shared_ptr<xolotlCore::ReactionNetwork> getSimpleReactionNetwork(const int maxClusterSize = 10,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry = std::make_shared<xolotlPerf::DummyHandlerRegistry>());

} /* end namespace testUtils */
#endif
