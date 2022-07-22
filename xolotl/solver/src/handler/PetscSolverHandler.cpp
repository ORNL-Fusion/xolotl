#include <xolotl/solver/handler/PetscSolverHandler.h>

namespace xolotl
{
namespace solver
{
namespace handler
{
PetscSolverHandler::PetscSolverHandler(
	NetworkType& _network, const options::IOptions& options) :
	SolverHandler(_network, options),
	fluxTimer(perfHandler->getTimer("Flux")),
	partialDerivativeTimer(perfHandler->getTimer("Partial Derivatives")),
	fluxCounter(perfHandler->getEventCounter("Flux")),
	partialDerivativeCounter(
		perfHandler->getEventCounter("Partial Derivatives"))
{
}

std::vector<PetscInt>
PetscSolverHandler::ConvertToPetscSparseFillMap(
	size_t dof, const core::network::IReactionNetwork::SparseFillMap& fillMap)
{
	// Determine number of non-zeros
	uint64_t nNonZeros = 0;
	for (auto const& currMapItem : fillMap) {
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

std::array<std::vector<PetscInt>, 2>
PetscSolverHandler::convertToCoordinateListPair(std::size_t dof,
	const core::network::IReactionNetwork::SparseFillMap& fillMap)
{
	auto nNonZeros =
		std::accumulate(fillMap.begin(), fillMap.end(), std::uint64_t{0},
			[](std::uint64_t r, auto&& kvp) { return r + kvp.second.size(); });

	std::array<std::vector<PetscInt>, 2> ret;
	ret[0].reserve(nNonZeros);
	ret[1].reserve(nNonZeros);

	for (auto i = 0; i < dof; ++i) {
		auto rowIter = fillMap.find(i);
		if (rowIter == fillMap.end()) {
			continue;
		}
		for (auto j : rowIter->second) {
			ret[0].push_back(i);
			ret[1].push_back(j);
		}
	}

	return ret;
}

std::vector<core::RowColPair>
PetscSolverHandler::convertToRowColPairList(std::size_t dof,
	const core::network::IReactionNetwork::SparseFillMap& fillMap)
{
	auto nNonZeros =
		std::accumulate(fillMap.begin(), fillMap.end(), std::uint64_t{0},
			[](std::uint64_t r, auto&& kvp) { return r + kvp.second.size(); });

    std::vector<core::RowColPair> ret;
	ret.reserve(nNonZeros);

	for (auto i = 0; i < dof; ++i) {
		auto rowIter = fillMap.find(i);
		if (rowIter == fillMap.end()) {
			continue;
		}
		for (auto j : rowIter->second) {
			ret.push_back({i,j});
		}
	}

	return ret;
}
} /* end namespace handler */
} /* end namespace solver */
} /* end namespace xolotl */
