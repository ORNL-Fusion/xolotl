#pragma once

#include <array>
#include <string>
#include <vector>

#include <boost/rational.hpp>

#include <xolotl/config.h>
#include <xolotl/options/IOptions.h>
#include <xolotl/util/Filesystem.h>

namespace xolotl
{
namespace interface
{
class NetworkHandlerClassGenerator
{
public:
	NetworkHandlerClassGenerator() = default;

	explicit NetworkHandlerClassGenerator(const fs::path& rnFile,
		const std::string& material, const std::string& userMaterial);

	void
	generate();

	void
	generateClusterGenerator();

	void
	generateReactions();

	void
	generateTraits();

	void
	generateNetwork();

	void
	generateNetworkHandler();

	void
	generateMaterialHandler();

	void
	generateBuild();

	const fs::path&
	getLibraryFileName() const
	{
		return _networkLibFile;
	}

private:
	void
	readNetworkFile();

	void
	generateReactionHeader();

	void
	generateReactionImpl();

	void
	generateClusterGeneratorHeader();

	void
	generateClusterGeneratorImpl();

	void
	generateNetworkHeader();

	void
	generateNetworkImpl();

	void
	writeCMakeLists();

	bool
	needToGenerate() const;

private:
	fs::path _rnFile;
	std::string _material;
	std::string _userMaterial;
	fs::path _genDir;
	fs::path _buildDir;
	std::string _baseName;
	std::string _speciesList;
	std::string _prodReaction;
	std::string _dissReaction;
	std::string _sinkReaction;
	std::string _clusterGenerator;
	std::string _reactionGenerator;
	std::string _reactionNetwork;
	std::string _networkHandler;
	std::string _materialHandler;
	fs::path _networkLibFile;
	fs::path _execFile;

	struct NetworkData
	{
		std::string label;
		std::string desc;
		double latticeParameter;
		double atomicVolume;
		std::string material;
	} _networkData;

	struct SpeciesData
	{
		enum class Type
		{
			impurity,
			vacancy,
			interstitial
		};
		std::string name;
		std::string label;
		Type type;
		AmountType maxPure{bad<AmountType>};

		struct Interaction
		{
			std::string species;
			boost::rational<AmountType> ratio;
		};
		std::vector<Interaction> interact;
	};
	std::vector<SpeciesData> _speciesData;
	std::unordered_map<std::string, const SpeciesData*> _speciesLabelMap;

	struct ClusterType
	{
		std::string label{};
		const SpeciesData* species{nullptr};
		std::array<std::string, 2> constituentLabels{};
		std::array<const SpeciesData*, 2> constituents{nullptr, nullptr};
	};

	struct ClusterData
	{
		ClusterType type;
		std::array<AmountType, 2> size;
		std::string radiusExpr;
		double migrationEnergy;
		double diffusionFactor;
	};
	std::vector<ClusterData> _clusterData;
	struct ClusterGroup
	{
		ClusterType type;
		std::vector<const ClusterData*> cases;
		const ClusterData* general{nullptr};
	};
	std::unordered_map<std::string, ClusterGroup> _clusterGroups;

	struct ReactionParts
	{
		std::vector<std::string> reactantLabels;
		std::vector<ClusterType> reactants;
		std::vector<std::string> productLabels;
		std::vector<ClusterType> products;
	};
	struct ReactionData
	{
		ReactionParts parts;
		std::string type;
		AmountType size;
		std::string bindingExpr;
	};
	std::vector<ReactionData> _reactionData;
	struct ReactionGroup
	{
		ReactionParts parts;
		std::vector<const ReactionData*> cases;
		const ReactionData* general{nullptr};
	};
	std::unordered_map<std::string, ReactionGroup> _productionReactionGroups;
	std::unordered_map<std::string, ReactionGroup> _dissociationReactionGroups;
};

void
generateNetworkHandler(options::IOptions* options);
} // namespace interface
} // namespace xolotl
