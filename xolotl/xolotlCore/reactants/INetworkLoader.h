#ifndef INETWORKLOADER_H
#define INETWORKLOADER_H

#include "IReactionNetwork.h"
#include <IOptions.h>

namespace xolotlCore {

/**
 *  This class is the interface for the netowrk loader.
 */
class INetworkLoader {

public:

	/**
	 * The destructor.
	 */
	virtual ~INetworkLoader() {
	}

	/**
	 * This operation specifies the inputstream from which cluster data should
	 * be loaded.
	 *
	 * @param stream The inputstream from which the cluster data should be
	 * loaded
	 */
	virtual void setInputstream(const std::shared_ptr<std::istream> stream) = 0;

	/**
	 * This operation will load the reaction network from the input file in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 *
	 * @return network The reaction network
	 */
	virtual std::shared_ptr<IReactionNetwork> load() = 0;

	/**
	 * This operation will generate the reaction network from options.
	 * The network will be empty if it can not be loaded.
	 *
	 * @param options The command line options
	 * @return network The reaction network
	 */
	virtual std::shared_ptr<IReactionNetwork> generate(IOptions &options) = 0;

	/**
	 * This operation will set the name of the file where to take the network from.
	 *
	 * @param name The name of the file
	 */
	virtual void setFilename(const std::string& name) = 0;

	/**
	 * This operation will get the name of the file where to take the network from.
	 *
	 * @return The name of the file
	 */
	virtual std::string getFilename() const = 0;

	/**
	 * This operation will set the reactions to dummy reactions.
	 */
	virtual void setDummyReactions() = 0;

};

}

#endif
