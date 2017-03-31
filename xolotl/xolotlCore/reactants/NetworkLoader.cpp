// Includes
#include "NetworkLoader.h"

// Namespaces
using namespace xolotlCore;

NetworkLoader::NetworkLoader() : networkStream(nullptr), handlerRegistry(nullptr),
		fileName(""),
		dummyReactions(false) {}

NetworkLoader::NetworkLoader(std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		networkStream(nullptr), handlerRegistry(registry), fileName(""),
		dummyReactions(false) {}

NetworkLoader::NetworkLoader(const std::shared_ptr<std::istream> stream,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> registry) :
		networkStream(stream), handlerRegistry(registry), fileName(""),
		dummyReactions(false) {}

void NetworkLoader::setInputstream(
		const std::shared_ptr<std::istream> stream) {
	networkStream = stream;

	return;
}

void NetworkLoader::setFilename (const std::string& name) {
	fileName = name;

	return;
}

std::string NetworkLoader::getFilename () const {
	return fileName;
}

void NetworkLoader::setDummyReactions () {
	dummyReactions = true;

	return;
}
