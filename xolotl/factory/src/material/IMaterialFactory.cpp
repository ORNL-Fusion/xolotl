#include <xolotl/factory/material/IMaterialFactory.h>
#include <xolotl/factory/material/W100MaterialFactory.h>
#include <xolotl/factory/material/W110MaterialFactory.h>
#include <xolotl/factory/material/W111MaterialFactory.h>
#include <xolotl/factory/material/W211MaterialFactory.h>
#include <xolotl/factory/material/FuelMaterialFactory.h>
#include <xolotl/factory/material/TRIDYNMaterialFactory.h>
#include <xolotl/factory/material/PulsedMaterialFactory.h>
#include <xolotl/factory/material/AlloyMaterialFactory.h>
#include <xolotl/factory/material/FeMaterialFactory.h>

namespace xolotl {
namespace factory {
namespace material {

static std::shared_ptr<IMaterialFactory> theMaterialFactory;

std::shared_ptr<IMaterialFactory> IMaterialFactory::createMaterialFactory(
		const options::Options &opts) {
	// Get the material type
	const auto &materialType = opts.getMaterial();
	// W100 case
	if (materialType == "W100")
		theMaterialFactory = std::make_shared<W100MaterialFactory>(opts);
	// W110 case
	else if (materialType == "W110")
		theMaterialFactory = std::make_shared<W110MaterialFactory>(opts);
	// W111 case
	else if (materialType == "W111")
		theMaterialFactory = std::make_shared<W111MaterialFactory>(opts);
	// W211 case
	else if (materialType == "W211")
		theMaterialFactory = std::make_shared<W211MaterialFactory>(opts);
	// Fuel case
	else if (materialType == "Fuel")
		theMaterialFactory = std::make_shared<FuelMaterialFactory>(opts);
	// TRIDYN case
	else if (materialType == "TRIDYN")
		theMaterialFactory = std::make_shared<TRIDYNMaterialFactory>(opts);
	// Pulsed case
	else if (materialType == "Pulsed")
		theMaterialFactory = std::make_shared<PulsedMaterialFactory>(opts);
	// Alloy case
	else if (materialType == "800H")
		theMaterialFactory = std::make_shared<AlloyMaterialFactory>(opts);
	// Fe case
	else if (materialType == "Fe")
		theMaterialFactory = std::make_shared<FeMaterialFactory>(opts);
	// The type is not supported
	else {
		throw std::string(
				"\nThe material type is not known: \"" + materialType + "\"");
	}

	return theMaterialFactory;
}

} // end namespace material
} // end namespace factory
} // end namespace xolotl
