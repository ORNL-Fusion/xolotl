#ifndef SIMPLEREACTIONNETWORK_H_
#define SIMPLEREACTIONNETWORK_H_

#include <Reactant.h>

namespace testUtils {

/**
 * This class creates a simple reaction network used for testing. It contains
 * 10 each of He, V and I clusters and twenty HeV clusters.
 *
 * It does not register itself as the ReactionNetwork for its clusters because
 * of limitations with shared_ptrs and "this." So, the
 * TestUtils::getSimpleReactionNetwork() operation should always be called to
 * insure that it is properly initialized.
 */
class SimpleReactionNetwork : public xolotlCore::ReactionNetwork {

public:
	//! Constructor
	SimpleReactionNetwork();

	//! Destructor
	virtual ~SimpleReactionNetwork();
};

/**
 * This operation creates a SimpleReactionNetwork and makes sure that it is
 * properly registered with the clusters it contains. This operation should
 * always be called instead of constructing a SimpleReactionNetwork manually.
 * @return The reaction network.
 */
std::shared_ptr<xolotlCore::ReactionNetwork> getSimpleReactionNetwork();

} /* end namespace testUtils */
#endif
