#include "xolotlSolver/solverhandler/PetscSolverHandler.h"

namespace xolotlSolver {

std::vector<PetscInt> PetscSolverHandler::ConvertToPetscSparseFillMap(
		size_t dof,
		const xolotlCore::experimental::IReactionNetwork::SparseFillMap &fillMap) {

	// Determine number of non-zeros
	uint64_t nNonZeros = 0;
	for (auto const &currMapItem : fillMap) {
		nNonZeros += currMapItem.second.size();
	}

	// Allocate a 1D vector of the correct size.
	std::vector<PetscInt> ret(nNonZeros + dof + 1);

	// For each row, determine the starting index value of its items,
	// and note the column indices of its non-zeros.
	size_t currIdx = dof + 1;
	for (auto i = 0; i < dof; ++i) {
		// Set starting index of current row's items.
		ret[i] = currIdx;

		// Fill out current row's items with column ids of non-zeros.
		auto rowIter = fillMap.find(i);
		if (rowIter != fillMap.end()) {
			// We have some non-zeros for current row.
			for (auto colId : rowIter->second) {
				ret[currIdx] = colId;
				++currIdx;
			}
		}
	}

	// Set the end of the array.
	ret[dof] = (nNonZeros + dof + 1);

	return ret;
}

} // nmaespace xolotlSolver
