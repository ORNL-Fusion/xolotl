#include <IMaterialFactory.h>
#include <W100MaterialFactory.h>
#include <W110MaterialFactory.h>
#include <W111MaterialFactory.h>
#include <W211MaterialFactory.h>
#include <FuelMaterialFactory.h>
#include <TRIDYNMaterialFactory.h>
#include <PulsedMaterialFactory.h>
#include <AlloyMaterialFactory.h>
#include <FeMaterialFactory.h>
#include <UZrMaterialFactory.h>

namespace xolotlFactory {

static std::shared_ptr<IMaterialFactory> theMaterialFactory;

std::shared_ptr<IMaterialFactory> IMaterialFactory::createMaterialFactory(
		const xolotlCore::Options &options) {
	// Get the material type
	const auto &materialType = options.getMaterial();
	// W100 case
	if (materialType == "W100")
		theMaterialFactory = std::make_shared<W100MaterialFactory>(options);
	// W110 case
	else if (materialType == "W110")
		theMaterialFactory = std::make_shared<W110MaterialFactory>(options);
	// W111 case
	else if (materialType == "W111")
		theMaterialFactory = std::make_shared<W111MaterialFactory>(options);
	// W211 case
	else if (materialType == "W211")
		theMaterialFactory = std::make_shared<W211MaterialFactory>(options);
	// Fuel case
	else if (materialType == "Fuel")
		theMaterialFactory = std::make_shared<FuelMaterialFactory>(options);
	// TRIDYN case
	else if (materialType == "TRIDYN")
		theMaterialFactory = std::make_shared<TRIDYNMaterialFactory>(options);
	// Pulsed case
	else if (materialType == "Pulsed")
		theMaterialFactory = std::make_shared<PulsedMaterialFactory>(options);
	// Alloy case
	else if (materialType == "800H")
		theMaterialFactory = std::make_shared<AlloyMaterialFactory>(options);
	// Fe case
	else if (materialType == "Fe")
		theMaterialFactory = std::make_shared<FeMaterialFactory>(options);
	// UZr case
	else if (materialType == "UZr")
		theMaterialFactory = std::make_shared<UZrMaterialFactory>(options);
	// The type is not supported
	else {
		throw std::string(
				"\nThe material type is not known: \"" + materialType + "\"");
	}

	return theMaterialFactory;
}

} // end namespace xolotlFactory
