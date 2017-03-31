#ifndef IREACTIONHANDLERFACTORY_H
#define IREACTIONHANDLERFACTORY_H

#include <Options.h>
#include <INetworkLoader.h>
#include <IReactionNetwork.h>

namespace xolotlFactory {

/**
 * Realizations of this interface are responsible for handling the flux and the advection.
 * they are both dependent on the type of material under study.
 */
class IReactionHandlerFactory {
public:

	/**
	 * The destructor
	 */
	~IReactionHandlerFactory() {}

	/**
	 * Initialize the reaction network.
	 *
	 * @param options The options.
	 * @param registry The performance registry.
	 */
	virtual void initializeReactionNetwork(xolotlCore::Options &options,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) = 0;

	/**
	 * Return the network loader.
	 *
	 * @return The network loader.
	 */
	virtual std::shared_ptr<xolotlCore::INetworkLoader> getNetworkLoaderHandler() const = 0;

	/**
	 * Return the network.
	 *
	 * @return The network.
	 */
	virtual std::shared_ptr<xolotlCore::IReactionNetwork> getNetworkHandler() const = 0;

	/**
	 * Function that create the wanted reaction handler factory depending on the given type.
	 *
	 * @param problemType The type of wanted problem (PSI or NE).
	 * @return The reaction factory.
	 */
	static std::shared_ptr<IReactionHandlerFactory> createNetworkFactory(const std::string& problemType);

};

} // end namespace xolotlFactory

#endif // IREACTIONHANDLERFACTORY_H
