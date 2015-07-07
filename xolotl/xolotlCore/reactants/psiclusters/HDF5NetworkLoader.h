#ifndef HDF5NETWORKLOADER_H_
#define HDF5NETWORKLOADER_H_

//Includes
#include <PSIClusterNetworkLoader.h>

namespace xolotlCore {

/**
 * This class overwrites the load() methods of PSIClusterNetworkLoader
 * for HDF5 files.
 */
class HDF5NetworkLoader: public PSIClusterNetworkLoader {
private:

	/**
	 * Name of the file to load the network from.
	 */
	std::string fileName;

	/**
	 * Private nullary constructor.
	 */
	HDF5NetworkLoader() {}

public:

	/**
	 * The default constructor. The setInputstream() operation must be called
	 * if this constructor is used.
	 */
	HDF5NetworkLoader(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) {
		handlerRegistry = registry;
	}

	/**
	 * The destructor.
	 */
	virtual ~HDF5NetworkLoader() {
	}

	/**
	 * This operation will load the reaction network from the HDF5 file in
	 * the format specified previously. The network will be empty if it can not
	 * be loaded.
	 *
	 * @return The reaction network.
	 */
	std::shared_ptr<PSIClusterReactionNetwork> load();

	/**
	 * This operation will set the name of the file where to take the network from.
	 *
	 * @param name The name of the file
	 */
	void setFilename (const std::string& name);

	/**
	 * This operation will get the name of the file where to take the network from.
	 *
	 * @return The name of the file
	 */
	std::string getFilename () const;
};

} /* namespace xolotlCore */

#endif /* HDF5NETWORKLOADER_H_ */
