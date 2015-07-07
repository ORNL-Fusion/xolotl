#include <IMaterialFactory.h>
#include <W100MaterialFactory.h>
#include <W110MaterialFactory.h>
#include <W111MaterialFactory.h>

namespace xolotlFactory {

static std::shared_ptr<IMaterialFactory> theMaterialFactory;

std::shared_ptr<IMaterialFactory> IMaterialFactory::createMaterialFactory(const std::string& materialType,
		int dimension) {
	// W100 case
	if (materialType == "W100") theMaterialFactory = std::make_shared<W100MaterialFactory>(dimension);
	// W110 case
	else if (materialType == "W110") theMaterialFactory = std::make_shared<W110MaterialFactory>(dimension);
	// W111 case
	else if (materialType == "W111") theMaterialFactory = std::make_shared<W111MaterialFactory>(dimension);
	// The type is not supported
	else {
		throw std::string(
				"\nThe material type is not known: \"" + materialType
				+ "\"");
	}
	return theMaterialFactory;
}

};  // end namespace xolotlFactory
