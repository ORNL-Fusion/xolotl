// Includes
#include <PetscSolverHandler.h>
#include <HDF5Utils.h>
#include <MathUtils.h>
#include <Constants.h>

namespace xolotlSolver {

// Definitions of some external variables (declared in PetscSolver1DHandler.h)
double lastTemperature = 0.0;
std::unordered_map<int, std::vector<int> > dFillMap;
std::shared_ptr<std::vector<xolotlCore::Reactant *>> allReactants;
std::vector<double> clusterPartials;

/**
 * A vector for holding the partial derivatives for one cluster in the order
 * that PETSc expects. It is sized in the createSolverContext() operation.
 *
 * The vector is used for every cluster and immediately reset to zero before
 * being used for the next. This allows the acquisition of the partial
 * derivatives to take up minimal memory and require no additional dynamic
 * allocations.
 */
std::vector<double> reactingPartialsForCluster;

void PetscSolverHandler::getDiagonalFill(PetscInt *diagFill,
		int diagFillSize) {

	// Degrees of freedom is the total number of clusters in the network
	const int dof = network->size();
	const int diagSize = dof * dof;

	// Declarations for the loop
	std::vector<int> connectivity;
	int connectivityLength, id, index;

	// Fill the diagonal block if the sizes match up
	if (diagFillSize == diagSize) {
		// Get the connectivity for each reactant
		for (int i = 0; i < dof; i++) {
			// Get the reactant and its connectivity
			auto reactant = allReactants->at(i);
			connectivity = reactant->getConnectivity();
			connectivityLength = connectivity.size();
			// Get the reactant id so that the connectivity can be lined up in
			// the proper column
			id = reactant->getId() - 1;
			// Create the vector that will be inserted into the dFill map
			std::vector<int> columnIds;
			// Add it to the diagonal fill block
			for (int j = 0; j < connectivityLength; j++) {
				// The id starts at j*connectivity length and is always offset
				// by the id, which denotes the exact column.
				index = id * connectivityLength + j;
				diagFill[index] = connectivity[j];
				// Add a column id if the connectivity is equal to 1.
				if (connectivity[j] == 1) {
					columnIds.push_back(j);
				}
			}
			// Update the map
			dFillMap[id] = columnIds;
		}
		// Debug output
//		std::cout << "Number of degrees of freedom = " << dof
//				<< std::endl;
//		printf("\n");
//		for (i = 0; i < dof; i++) {
//			for (j = 0; j < dof; j++) {
//				printf("%d ", dfill[i * dof + j]);
//			}
//			printf("\n");
//		}
//		printf("\n");
	} else {
		std::string err =
				"PetscSolverHandler Exception: Invalid diagonal block size!\n";
		throw std::string(err);
	}

	return;
}

} /* end namespace xolotlSolver */
