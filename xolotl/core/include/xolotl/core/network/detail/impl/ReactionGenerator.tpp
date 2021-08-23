#pragma once

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
template <typename TNetwork, typename TDerived>
ReactionGeneratorBase<TNetwork, TDerived>::ReactionGeneratorBase(
	const TNetwork& network) :
	_subpaving(network._subpaving),
	_clusterData(network._clusterData.h_view()),
	_clusterDataView(network._clusterData.d_view),
	_numDOFs(network.getDOF()),
	_enableReducedJacobian(network.getEnableReducedJacobian()),
	_clusterProdReactionCounts(
		"Production Reaction Counts", _clusterData.numClusters),
	_clusterDissReactionCounts(
		"Dissociation Reaction Counts", _clusterData.numClusters)
{
}

template <typename TNetwork, typename TDerived>
ReactionCollection<TNetwork>
ReactionGeneratorBase<TNetwork, TDerived>::generateReactions()
{
	auto numClusters = _clusterData.numClusters;
	auto diffusionFactor = _clusterData.diffusionFactor;
	auto generator = *(this->asDerived());
	using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2>>;
	auto range2d = Range2D({0, 0}, {numClusters, numClusters});
	Kokkos::parallel_for(
		range2d, KOKKOS_LAMBDA(IndexType i, IndexType j) {
			if (j < i) {
				return;
			}
			if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
				return;
			}
			generator(i, j, Count{});
		});
	Kokkos::fence();

	setupCrs();

	generator = *(this->asDerived());

	Kokkos::parallel_for(
		range2d, KOKKOS_LAMBDA(IndexType i, IndexType j) {
			if (j < i) {
				return;
			}
			if (diffusionFactor(i) == 0.0 && diffusionFactor(j) == 0.0) {
				return;
			}
			generator(i, j, Construct{});
		});
	Kokkos::fence();

	// TODO: Should this be done in the ReactionCollection constructor?
	//      - Constructing all reactions
	//      - Generating connectivity
	auto reactionCollection = this->asDerived()->getReactionCollection();
	reactionCollection.constructAll(_clusterDataView, _allClusterSets);

	Kokkos::fence();

	generateConnectivity(reactionCollection);

	return reactionCollection;
}

template <typename TNetwork, typename TDerived>
typename TNetwork::IndexType
ReactionGeneratorBase<TNetwork, TDerived>::getRowMapAndTotalReactionCount()
{
	_numProdReactions = Kokkos::get_crs_row_map_from_counts(
		_prodCrsRowMap, _clusterProdReactionCounts);
	_numDissReactions = Kokkos::get_crs_row_map_from_counts(
		_dissCrsRowMap, _clusterDissReactionCounts);

	_prodReactions = Kokkos::View<ProductionReactionType*>(
		"Production Reactions", _numProdReactions);
	_dissReactions = Kokkos::View<DissociationReactionType*>(
		"Dissociation Reactions", _numDissReactions);

	return _numProdReactions + _numDissReactions;
}

template <typename TNetwork, typename TDerived>
void
ReactionGeneratorBase<TNetwork, TDerived>::setupCrsClusterSetSubView()
{
	_prodCrsClusterSets = getClusterSetSubView(
		std::make_pair(static_cast<IndexType>(0), _numProdReactions));

	_dissCrsClusterSets = getClusterSetSubView(std::make_pair(
		_numProdReactions, _numProdReactions + _numDissReactions));
}

template <typename TNetwork, typename TDerived>
void
ReactionGeneratorBase<TNetwork, TDerived>::setupCrs()
{
	auto numTotalReactions =
		this->asDerived()->getRowMapAndTotalReactionCount();
	_allClusterSets = ClusterSetView("Cluster Sets", numTotalReactions);
	this->asDerived()->setupCrsClusterSetSubView();
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReactionGeneratorBase<TNetwork, TDerived>::addProductionReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!_clusterData.enableStdReaction())
		return;

	Kokkos::atomic_increment(&_clusterProdReactionCounts(clusterSet.cluster0));
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReactionGeneratorBase<TNetwork, TDerived>::addProductionReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!_clusterData.enableStdReaction())
		return;

	auto id = _prodCrsRowMap(clusterSet.cluster0);
	for (; !Kokkos::atomic_compare_exchange_strong(
			 &_prodCrsClusterSets(id).cluster0, NetworkType::invalidIndex(),
			 clusterSet.cluster0);
		 ++id) { }
	_prodCrsClusterSets(id) = clusterSet;
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReactionGeneratorBase<TNetwork, TDerived>::addDissociationReaction(
	Count, const ClusterSet& clusterSet) const
{
	if (!_clusterData.enableStdReaction())
		return;

	Kokkos::atomic_increment(&_clusterDissReactionCounts(clusterSet.cluster1));
}

template <typename TNetwork, typename TDerived>
KOKKOS_INLINE_FUNCTION
void
ReactionGeneratorBase<TNetwork, TDerived>::addDissociationReaction(
	Construct, const ClusterSet& clusterSet) const
{
	if (!_clusterData.enableStdReaction())
		return;

	auto id = _dissCrsRowMap(clusterSet.cluster1);
	for (; !Kokkos::atomic_compare_exchange_strong(
			 &_dissCrsClusterSets(id).cluster1, NetworkType::invalidIndex(),
			 clusterSet.cluster1);
		 ++id) { }
	_dissCrsClusterSets(id) = clusterSet;
}

template <typename TNetwork, typename TDerived>
void
ReactionGeneratorBase<TNetwork, TDerived>::generateConnectivity(
	ReactionCollection<NetworkType>& reactionCollection)
{
	using RowMap = typename Connectivity::row_map_type;
	using Entries = typename Connectivity::entries_type;

	Connectivity tmpConn;
	// Count connectivity entries
	// NOTE: We're using row_map for counts because
	//      Reaction::contributeConnectivity expects the connectivity CRS
	tmpConn.row_map = RowMap(
		Kokkos::ViewAllocateWithoutInitializing("tmp counts"), this->_numDOFs);
	// Even if there is no reaction each dof should connect with itself (for
	// PETSc)
	Kokkos::parallel_for(
		this->_numDOFs,
		KOKKOS_LAMBDA(const IndexType i) { tmpConn.row_map(i) = 1; });
	if (this->_enableReducedJacobian) {
		reactionCollection.forEach(DEVICE_LAMBDA(auto&& reaction) {
			reaction.contributeReducedConnectivity(tmpConn);
		});
	}
	else {
		reactionCollection.forEach(DEVICE_LAMBDA(
			auto&& reaction) { reaction.contributeConnectivity(tmpConn); });
	}

	Kokkos::fence();
	// Get row map
	auto counts = tmpConn.row_map;
	auto nEntries =
		Kokkos::get_crs_row_map_from_counts(tmpConn.row_map, counts);
	// Reset counts view
	counts = RowMap();
	// Initialize entries to invalid
	tmpConn.entries =
		Entries(Kokkos::ViewAllocateWithoutInitializing("connectivity entries"),
			nEntries);
	Kokkos::parallel_for(
		nEntries, KOKKOS_LAMBDA(IndexType i) {
			tmpConn.entries(i) = NetworkType::invalidIndex();
		});
	// Even if there is no reaction each dof should connect with itself (for
	// PETSc)
	Kokkos::parallel_for(
		this->_numDOFs, KOKKOS_LAMBDA(const IndexType i) {
			auto id = tmpConn.row_map(i);
			for (; !Kokkos::atomic_compare_exchange_strong(
					 &tmpConn.entries(id), NetworkType::invalidIndex(), i);
				 ++id) {
				if (tmpConn.entries(id) == i) {
					break;
				}
			}
		});
	// Fill entries (column ids)
	if (this->_enableReducedJacobian) {
		reactionCollection.forEach(DEVICE_LAMBDA(auto&& reaction) {
			reaction.contributeReducedConnectivity(tmpConn);
		});
	}
	else {
		reactionCollection.forEach(DEVICE_LAMBDA(
			auto&& reaction) { reaction.contributeConnectivity(tmpConn); });
	}
	Kokkos::fence();

	// Shrink to fit
	Connectivity connectivity;
	Kokkos::count_and_fill_crs(
		connectivity, this->_numDOFs,
		KOKKOS_LAMBDA(IndexType i, IndexType * fill) {
			IndexType ret = 0;
			if (fill == nullptr) {
				auto jStart = tmpConn.row_map(i);
				auto jEnd = tmpConn.row_map(i + 1);
				ret = jEnd - jStart;
				for (IndexType j = jStart; j < jEnd; ++j) {
					if (tmpConn.entries(j) == NetworkType::invalidIndex()) {
						ret = j - jStart;
						break;
					}
				}
			}
			else {
				auto tmpStart = tmpConn.row_map(i);
				for (IndexType j = tmpStart; j < tmpConn.row_map(i + 1); ++j) {
					auto entry = tmpConn.entries(j);
					if (entry == NetworkType::invalidIndex()) {
						break;
					}
					fill[j - tmpStart] = entry;
				}
			}
			return ret;
		});
	nEntries = connectivity.entries.extent(0);

	_connectivity = connectivity;
	reactionCollection.setConnectivity(_connectivity);
}
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
