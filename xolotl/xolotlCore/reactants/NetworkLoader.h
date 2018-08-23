#ifndef NETWORKLOADER_H
#define NETWORKLOADER_H

#include "INetworkLoader.h"
#include <IHandlerRegistry.h>

namespace xolotlCore {

/**
 *  This class is the simplest implementation of the network loader.
 */
class NetworkLoader: public INetworkLoader {

protected:

	/**
	 * The istream from which the network of clusters will be read.
	 */
	std::shared_ptr<std::istream> networkStream;

	/**
	 * The performance handler registry used to measure runtime performance
	 * during loads.
	 */
	std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry;

	/**
	 * Name of the file to load the network from.
	 */
	std::string fileName;

	/**
	 * If we want dummy reactions or not
	 */
	bool dummyReactions;

	/**
	 * Private nullary constructor.
	 */
	NetworkLoader();

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 *
	 * @param registry The performance handler registry
	 */
	NetworkLoader(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * An alternative constructor provided for convenience.
	 *
	 * @param stream The inputstream from which the cluster data should be
	 * loaded
	 * @param registry The performance handler registry
	 */
	NetworkLoader(const std::shared_ptr<std::istream> stream,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry);

	/**
	 * The destructor.
	 */
	virtual ~NetworkLoader() {
	}

	/**
	 * This operation specifies the inputstream from which cluster data should
	 * be loaded.
	 *
	 * @param stream The inputstream from which the cluster data should be
	 * loaded
	 */
	void setInputstream(const std::shared_ptr<std::istream> stream);

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
	void setFilename(const std::string& name);

	/**
	 * This operation will get the name of the file where to take the network from.
	 *
	 * @return The name of the file
	 */
	std::string getFilename() const;

	/**
	 * This operation will set the reactions to dummy reactions.
	 */
	void setDummyReactions();

};

}

#endif
