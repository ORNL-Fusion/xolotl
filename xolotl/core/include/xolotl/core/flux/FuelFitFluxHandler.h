#pragma once

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/NEReactionNetwork.h>
#include <xolotl/util/MathUtils.h>
#include <xolotl/util/Tokenizer.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the FluxHandler interface to calculate the incident
 * xenon flux for nuclear fuel.
 */
class FuelFitFluxHandler : public FluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
	{
		// Constant flux
		return 1.0;
	}

	/**
	 * The fission yield for Xe.
	 */
	double xeYield;

	/**
	 * The fission yield for defects.
	 */
	double defectYield;

	// The temperature
	double temperature;

	// The list of pure vacancy indices
	std::vector<IdType> pureDefectIds;
	std::vector<double> pureDefectFactors;

public:
	/**
	 * The constructor
	 */
	FuelFitFluxHandler(const options::IOptions& options) :
		FluxHandler(options),
		xeYield(options.getFissionYield()),
		defectYield(1.0e4),
		temperature(options.getTempParam(0))
	{
	}

	/**
	 * The Destructor
	 */
	~FuelFitFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		// Set the grid
		xGrid = grid;

		// Skip if the flux amplitude is 0.0 and we are not using a time profile
		if (util::equal(fluxAmplitude, 0.0) && !useTimeProfile)
			return;

		// Set the flux index corresponding the the single xenon cluster here
		using NetworkType = network::NEReactionNetwork;
		auto& neNetwork = dynamic_cast<NetworkType&>(network);
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::V] = 1;
		auto cluster = neNetwork.findCluster(comp, plsm::HostMemSpace{});
		// Check that the helium cluster is present in the network
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::string(
				"\nThe single vacancy cluster is not present in the network, "
				"cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		pureDefectIds.push_back(cluster.getId());

		comp[NetworkType::Species::V] = 2;
		cluster = neNetwork.findCluster(comp, plsm::HostMemSpace{});
		// Check that the helium cluster is present in the network
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::string("\nThe di-vacancy cluster is not "
							  "present in the network, "
							  "cannot use the flux option!");
		}
		pureDefectIds.push_back(cluster.getId());

		comp[NetworkType::Species::V] = 0;
		comp[NetworkType::Species::I] = 1;
		cluster = neNetwork.findCluster(comp, plsm::HostMemSpace{});
		// Check that the helium cluster is present in the network
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::string("\nThe single interstitial cluster is not "
							  "present in the network, "
							  "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		comp[NetworkType::Species::I] = 0;
		comp[NetworkType::Species::Xe] = 1;
		cluster = neNetwork.findCluster(comp, plsm::HostMemSpace{});
		// Check that the helium cluster is present in the network
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::string("\nThe single xenon cluster is not "
							  "present in the network, "
							  "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());

		// Read the information from the text file
		pureDefectFactors = std::vector<double>(2, 0.0);
		constexpr double k_B = ::xolotl::core::kBoltzmann;
		std::ifstream reactionFile;
		reactionFile.open("reactionRates.txt");
		// Get the line
		std::string line;
		getline(reactionFile, line);
		// Read the first line
		std::vector<double> tokens;
		util::Tokenizer<double>{line}(tokens);
		// And start looping on the lines
		while (tokens.size() > 0) {
			// Find the Id of the cluster
			NetworkType::Composition comp = NetworkType::Composition::zero();
			comp[NetworkType::Species::Xe] = static_cast<IdType>(tokens[0]);
			comp[NetworkType::Species::V] = static_cast<IdType>(tokens[1]);
			comp[NetworkType::Species::I] = static_cast<IdType>(tokens[2]);

			if (comp[NetworkType::Species::Xe] > 0 or
				comp[NetworkType::Species::I] > 0) {
				getline(reactionFile, line);
				if (line == "Reactions")
					break;
				tokens = util::Tokenizer<double>{line}();
				continue;
			}

			// Get the largetG0
			double largestG0 = tokens[3] - k_B * temperature * tokens[4];
			//			for (auto i = 3; i < tokens.size(); i += 4) {
			//				auto g0 = tokens[i] - k_B * temperature * tokens[i +
			// 1]; 				if (g0 < largestG0) largestG0 = g0;
			//			}
			// Loop on the linked clusters
			for (auto i = 3; i < tokens.size(); i += 4) {
				// Get its properties
				auto g0 = tokens[i] - k_B * temperature * tokens[i + 1];
				pureDefectFactors[comp[NetworkType::Species::V] - 1] +=
					comp[NetworkType::Species::V] *
					std::exp((largestG0 - g0) / (k_B * temperature));
			}

			getline(reactionFile, line);
			if (line == "Reactions")
				break;

			tokens = util::Tokenizer<double>{line}();
		}

		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	computeIncidentFlux(double currentTime, double* concOffset,
		double* updatedConcOffset, int xi, int surfacePos)
	{
		// Skip if no index was set
		if (fluxIndices.size() == 0)
			return;

		// Compute the available site fraction
		double availVFraction = 1.0; // nm-3
		for (auto i = 0; i < pureDefectFactors.size(); i++) {
			availVFraction -=
				concOffset[pureDefectIds[i]] * pureDefectFactors[i];
		}

		updatedConcOffset[fluxIndices[0]] +=
			defectYield * fluxAmplitude * availVFraction; // V
		updatedConcOffset[fluxIndices[1]] +=
			defectYield * fluxAmplitude * availVFraction; // I
		updatedConcOffset[fluxIndices[2]] += xeYield * fluxAmplitude; // Xe

		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	virtual void
	setFissionYield(double yield)
	{
		xeYield = yield;
	}
};
// end class FuelFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl
