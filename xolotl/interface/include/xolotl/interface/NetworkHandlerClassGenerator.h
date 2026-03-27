#pragma once

#include <array>
#include <string>
#include <vector>

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
		std::array<AmountType, 2> bounds;
	};
	std::vector<SpeciesData> _speciesData;
	std::unordered_map<std::string, const SpeciesData*> _speciesLabelMap;

	struct ClusterData
	{
		std::string type;
		std::array<AmountType, 2> size;
		std::string radiusExpr;
		double migrationEnergy;
		double diffusionFactor;
	};
	std::vector<ClusterData> _clusterData;
	std::vector<std::vector<const ClusterData*>> _clusterGroups;

	struct ReactionData
	{
		std::vector<std::string> reactantLabels;
		std::vector<const SpeciesData*> reactants;
		std::vector<std::string> productLabels;
		std::vector<const SpeciesData*> products;
		std::string type;
		std::string bindingExpr;
	};
	std::vector<ReactionData> _reactionData;
	struct ReactionGroup
	{
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
