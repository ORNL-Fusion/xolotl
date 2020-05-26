#ifndef IREACTIONHANDLERFACTORY_H
#define IREACTIONHANDLERFACTORY_H

#include <Options.h>
#include <experimental/IReactionNetwork.h>

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
    virtual ~IReactionHandlerFactory() {
	}

	/**
	 * Initialize the reaction network.
	 *
	 * @param options The options.
	 * @param registry The performance registry.
	 */
	virtual void initializeReactionNetwork(const xolotlCore::Options &options,
			std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) = 0;

	/**
	 * Return the network.
	 *
	 * @return The network.
	 */
	virtual xolotlCore::experimental::IReactionNetwork& getNetworkHandler() const = 0;

	/**
	 * Function that create the wanted reaction handler factory depending on the given type.
	 *
	 * @param problemType The type of wanted problem (PSI, NE, or Fe).
	 * @return The reaction factory.
	 */
	static std::shared_ptr<IReactionHandlerFactory> createNetworkFactory(
			const std::string& problemType);

    static void resetNetworkFactory();
};

} // end namespace xolotlFactory

#endif // IREACTIONHANDLERFACTORY_H
