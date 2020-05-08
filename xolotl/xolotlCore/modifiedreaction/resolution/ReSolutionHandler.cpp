// Includes
#include <ReSolutionHandler.h>

namespace xolotlCore {

void ReSolutionHandler::initialize(const IReactionNetwork& network,
		double electronicStoppingPower) {
	// The daughter class will do this part
	return;
}

void ReSolutionHandler::updateReSolutionRate(double rate) {

	return;
}

void ReSolutionHandler::setFissionYield(double yield) {
	fissionYield = yield;

	return;
}

void ReSolutionHandler::computeReSolution(const IReactionNetwork& network,
		double *concOffset, double *updatedConcOffset, int xi, int xs, int yj,
		int zk) {
	return;
}

int ReSolutionHandler::computePartialsForReSolution(
		const IReactionNetwork& network, double *val, int *indices, int xi,
		int xs, int yj, int zk) {
	return 0;
}

}/* end namespace xolotlCore */

