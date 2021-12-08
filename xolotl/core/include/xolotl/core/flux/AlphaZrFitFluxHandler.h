#ifndef ALPHAZRFITFLUXHANDLER_H
#define ALPHAZRFITFLUXHANDLER_H

#include <cmath>

#include <xolotl/core/flux/FluxHandler.h>
#include <xolotl/core/network/ZrReactionNetwork.h>

namespace xolotl
{
namespace core
{
namespace flux
{
/**
 * This class realizes the IFluxHandler interface to calculate the incident
 * fluxes for an alpha Zr material.
 */
class AlphaZrFitFluxHandler : public FluxHandler
{
private:
	/**
	 * \see FluxHandler.h
	 */
	double
	FitFunction(double x)
	{
		// Not actually used
		return 0.0;
	}

	// Define the range of possible cluster sizes, and their respective
	// fluxes in nm^-3 s^-1 Multiply by 1.0e+21 to get to cm^-3 s^-1 if
	// you want to check values
	float fluxI[70] = {2.365739684679306e-07, 2.1755974706776704e-08,
		6.6215695655160775e-09, 3.979857734281858e-09, 2.108954194096771e-09,
		1.7928725927606022e-09, 1.5213778723816244e-09, 5.500284627584542e-10,
		4.093585813506016e-10, 6.703070630459465e-10, 4.3671599902652566e-10,
		1.5136519821278934e-10, 2.0253029017681254e-10, 1.4309830185036925e-10,
		3.395551051218416e-10, 1.6644273710076962e-10, 9.872932486119687e-11,
		6.210012288223767e-11, 1.741549602054323e-10, 1.1421807857184855e-10,
		7.142901860019398e-11, 1.9311062910580072e-11, 1.389313676998744e-10,
		4.5282886983331086e-11, 2.8966594365870508e-11, 1.115499440160588e-10,
		0.0, 3.562735552804218e-11, 9.655531455290036e-12,
		3.562735552804218e-11, 9.655531455290036e-12, 0.0,
		2.597182407275075e-11, 2.597182407275075e-11, 0.0,
		2.597182407275075e-11, 0.0, 9.655531455290036e-12, 0.0, 0.0,
		9.655531455290036e-12, 0.0, 2.597182407275075e-11, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	float fluxV[70] = {1.7811390508204618e-07, 1.6930667065081818e-08,
		7.0735654428151215e-09, 3.3754591361812474e-09, 1.840000680680739e-09,
		1.6146127259999994e-09, 2.1627167789999997e-09, 1.015531922447914e-09,
		6.790095549999998e-10, 6.365080860000001e-10, 5.00133545e-10,
		6.908325919999998e-10, 6.111473429999999e-10, 3.894946049999999e-10,
		2.2615779300000004e-10, 4.01893729e-10, 1.45108826e-10,
		1.6787280799999993e-10, 2.0438660999999995e-10, 1.6719651500000004e-10,
		1.4427961800000002e-10, 1.61956032e-10, 9.429843900000001e-11,
		7.504127899999999e-11, 4.2669623e-11, 1.01114429e-10,
		1.5158115299999994e-10, 3.5527909e-11, 3.3041043e-11,
		1.2437397699999997e-10, 0, 7.1229639e-11, 6.856895199999998e-11,
		4.515648899999999e-11, 6.731392999999998e-11, 3.3041043e-11, 0,
		8.501352199999999e-11, 2.8885739999999996e-11, 0, 3.5527909e-11, 0,
		5.1972479000000005e-11, 0, 9.628580000000001e-12, 0,
		9.628580000000001e-12, 1.9257160000000002e-11, 1.9257160000000002e-11,
		3.3041043e-11, 9.628580000000001e-12, 1.9257160000000002e-11, 0, 0, 0,
		2.5899328999999996e-11, 2.5899328999999996e-11, 3.428447599999999e-11,
		9.628580000000001e-12, 0, 0, 0, 2.5899328999999996e-11, 0, 0,
		1.9257160000000002e-11, 0, 0, 0, 0};

	// Set the maximum cluster sizes to set a generation flux to
	int maxSizeI = 70;
	int maxSizeV = 70;

public:
	/**
	 * The constructor
	 */
	AlphaZrFitFluxHandler(const options::IOptions& options) :
		FluxHandler(options)
	{
	}

	/**
	 * The Destructor
	 */
	~AlphaZrFitFluxHandler()
	{
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	initializeFluxHandler(network::IReactionNetwork& network, int surfacePos,
		std::vector<double> grid)
	{
		// std::cout << "Flux-start \n";
		using NetworkType = network::ZrReactionNetwork;
		auto zrNetwork = dynamic_cast<NetworkType*>(&network);

		// Set the flux index corresponding the mobile interstitial clusters (n
		// < 10)
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::I] = 1;
		std::ostringstream oss;
		auto cluster = zrNetwork->findCluster(comp, plsm::onHost);
		for (int i = 1; i < (maxSizeI + 1); i++) {
			comp[NetworkType::Species::I] = i;
			cluster = zrNetwork->findCluster(comp, plsm::onHost);
			if (cluster.getId() == NetworkType::invalidIndex()) {
				throw std::runtime_error(
					"\nThe current interstitial cluster is not "
					"present in the network, "
					"cannot use the flux option!");
			}
			fluxIndices.push_back(cluster.getId());
			incidentFluxVec.push_back(std::vector<double>(1, fluxI[i - 1]));
		}

		// Set the flux index corresponding the mobile vacancy clusters (n < 10)
		comp[NetworkType::Species::I] = 0;
		for (int i = 1; i < (maxSizeV + 1); i++) {
			comp[NetworkType::Species::V] = i;
			cluster = zrNetwork->findCluster(comp, plsm::onHost);
			if (cluster.getId() == NetworkType::invalidIndex()) {
				throw std::runtime_error("\nThe current vacancy cluster is not "
										 "present in the network, "
										 "cannot use the flux option!");
			}
			fluxIndices.push_back(cluster.getId());
			incidentFluxVec.push_back(std::vector<double>(1, fluxV[i - 1]));
		}

		/*
		comp[NetworkType::Species::V] = 0;
		comp[NetworkType::Species::Basal] = 1;
		cluster = zrNetwork->findCluster(comp, plsm::onHost);
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe current basal cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(std::vector<double>(1, fluxV[0]));


		comp[NetworkType::Species::Basal] = 2;
		cluster = zrNetwork->findCluster(comp, plsm::onHost);
		if (cluster.getId() == NetworkType::invalidIndex()) {
			throw std::runtime_error("\nThe current basal cluster is not "
									 "present in the network, "
									 "cannot use the flux option!");
		}
		fluxIndices.push_back(cluster.getId());
		incidentFluxVec.push_back(std::vector<double>(1, fluxV[1]));
		 */
		// std::cout << "Flux-stop \n";
		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	computeIncidentFlux(
		double currentTime, double* updatedConcOffset, int xi, int surfacePos)
	{
		// Define only for a 0D case
		if (xGrid.size() == 0) {
			// std::cout << "Flux-start2 \n";

			for (int i = 0; i < fluxIndices.size(); i++) {
				updatedConcOffset[fluxIndices[i]] += incidentFluxVec[i][0];
			}
		}

		else {
			throw std::runtime_error(
				"\nThe alpha Zr problem is not defined for more than 0D!");
		}

		return;
	}

	/**
	 * \see IFluxHandler.h
	 */
	std::vector<std::pair<IdType, double>>
	getImplantedFlux(std::vector<AmountType> map)
	{
		std::vector<std::pair<IdType, double>> toReturn;
		// Loop on the map
		for (auto i = 0; i < map.size(); i++) {
			// Look for this value in fluxIndices
			for (auto j = 0; i < fluxIndices.size(); j++) {
				if (map[i] == fluxIndices[j]) {
					toReturn.push_back(
						std::make_pair(i, incidentFluxVec[j][0]));
					break;
				}
			}
		}
		return toReturn;
	}

	/**
	 * \see IFluxHandler.h
	 */
	void
	setImplantedFlux(std::vector<std::pair<IdType, double>> fluxVector)
	{
		fluxIndices.clear();
		incidentFluxVec.clear();

		// Loop on the flux vector
		for (auto pair : fluxVector) {
			fluxIndices.push_back(pair.first);
			incidentFluxVec.push_back(std::vector<double>(1, pair.second));
		}
		return;
	}
};
// end class AlphaZrFitFluxHandler

} // namespace flux
} // namespace core
} // namespace xolotl

#endif
