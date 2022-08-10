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

	std::vector<double> fluxI = {4.141478088722051e-07, 3.961706397610119e-08,
		1.2106952158361212e-08, 7.648435047264139e-09, 4.0244599086078205e-09,
		3.3920182333877126e-09, 2.973973662425018e-09, 1.2417969019237562e-09,
		8.851070516960898e-10, 1.3603816064935588e-09, 7.755269088868743e-10,
		2.8874033916240647e-10, 4.672779143981795e-10, 3.2236846376884343e-10,
		7.656476570683874e-10, 3.573797643046147e-10, 2.1021120447778106e-10,
		1.539555427405712e-10, 3.782044331935486e-10, 2.568277326339155e-10,
		1.585689463091986e-10, 5.5696756415877455e-11, 3.3550319027591227e-10,
		1.069608242891428e-10, 8.354513462381556e-11, 2.4778164052265936e-10,
		0.0, 7.911244608120112e-11, 2.7848378207938727e-11,
		7.911244608120112e-11, 2.7848378207938727e-11, 0.0,
		5.1264067873265366e-11, 5.1264067873265366e-11, 0.0,
		5.1264067873265366e-11, 0.0, 2.7848378207938727e-11, 0.0, 0.0,
		2.7848378207938727e-11, 0.0, 5.1264067873265366e-11, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

	std::vector<double> fluxV = {3.000166061596907e-07, 3.261617020215222e-08,
		1.125375463678962e-08, 6.367679929976148e-09, 3.6816615428727384e-09,
		2.9269034006169983e-09, 3.359947468090799e-09, 1.8195145694037134e-09,
		1.3135561927353992e-09, 1.0449271185930001e-09, 1.1127620938895e-09,
		1.474853409405799e-09, 1.1821355154481997e-09, 8.466205016950002e-10,
		4.427436436175002e-10, 8.012900146363999e-10, 3.0815704976440003e-10,
		3.241739333764998e-10, 4.196949035263e-10, 3.490487822005e-10,
		3.084457042192001e-10, 3.664849375693001e-10, 2.1650475950509992e-10,
		1.625140160880999e-10, 9.524908203469999e-11, 2.125413615481001e-10,
		3.2572815399419994e-10, 7.668920153209999e-11, 6.825371032619996e-11,
		2.292135057965999e-10, 0, 1.5371191870060005e-10,
		1.4494291185829992e-10, 1.036845732406e-10, 1.7151468025099992e-10,
		6.825371032619996e-11, 0, 1.6797488560979994e-10, 8.098611512549996e-11,
		0, 7.668920153209999e-11, 0, 9.972117528360005e-11, 0,
		2.699537170850001e-11, 0, 2.699537170850001e-11, 5.399074341700002e-11,
		5.399074341700002e-11, 6.825371032619996e-11, 2.699537170850001e-11,
		5.399074341700002e-11, 0, 0, 0, 4.9693829823599966e-11,
		4.9693829823599966e-11, 1.4494291185829992e-10, 2.699537170850001e-11,
		0, 0, 0, 4.9693829823599966e-11, 0, 0, 5.399074341700002e-11, 0, 0, 0,
		0};

	// Keep the maximum cluster sizes to set a generation flux to
	size_t maxSizeI = 0;
	size_t maxSizeV = 0;
	size_t maxSizeB = 0;
	double Qb = 0;

public:
	/**
	 * The constructor
	 */
	AlphaZrFitFluxHandler(const options::IOptions& options) :
		FluxHandler(options),
		maxSizeI(options.getMaxI()),
		maxSizeV(options.getMaxV()),
		maxSizeB(options.getMaxImpurity())
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
		using NetworkType = network::ZrReactionNetwork;
		auto zrNetwork = dynamic_cast<NetworkType*>(&network);

        // Set the fraction of large vacancy clusters (n > 19) that become faulted basal pyramids:
        if (maxSizeB > 18) Qb = 0.10; // Basal
        else Qb = 0; //No basal

		// Set the flux index corresponding the mobile interstitial clusters (n
		// < 10)
		NetworkType::Composition comp = NetworkType::Composition::zero();
		comp[NetworkType::Species::I] = 1;
		std::ostringstream oss;
		auto cluster = zrNetwork->findCluster(comp, plsm::HostMemSpace{});
		for (int i = 1; i <= std::min(maxSizeI, fluxI.size()); i++) {
			comp[NetworkType::Species::I] = i;
			cluster = zrNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			incidentFluxVec.push_back(std::vector<double>(1, fluxI[i - 1]));
		}

		// Set the flux index corresponding the mobile vacancy clusters (n < 10)
		comp[NetworkType::Species::I] = 0;
		for (int i = 1; i <= std::min(maxSizeV, fluxV.size()); i++) {
			comp[NetworkType::Species::V] = i;
			cluster = zrNetwork->findCluster(comp, plsm::HostMemSpace{});
			if (cluster.getId() == NetworkType::invalidIndex()) {
				continue;
			}
			fluxIndices.push_back(cluster.getId());
			if (i >= 19)
				incidentFluxVec.push_back(
					std::vector<double>(1, fluxV[i - 1] * (1 - Qb)));
			else
				incidentFluxVec.push_back(std::vector<double>(1, fluxV[i - 1]));
		}

		// Set the flux index corresponding to Basal
		comp[NetworkType::Species::V] = 0;
		for (int i = 1; i <= std::min(maxSizeB, fluxV.size()); i++) {
			comp[NetworkType::Species::Basal] = i;
			cluster = zrNetwork->findCluster(comp, plsm::HostMemSpace{});
            if (cluster.getId() == NetworkType::invalidIndex()) {
                continue;
            }
            fluxIndices.push_back(cluster.getId());
            if (i >= 19)
                incidentFluxVec.push_back(std::vector<double>(1, fluxV[i - 1] * (Qb)));
            else
                incidentFluxVec.push_back(std::vector<double>(1, 0));
        }

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

            double cascadeEfficiency = (0.495*(1-tanh(0.00040527088*(currentTime/100-5000)))+0.025); //THIS ONE without basal now apparently?
            
            for (int i = 0; i < fluxIndices.size(); i++) {
                updatedConcOffset[fluxIndices[i]] += incidentFluxVec[i][0] * cascadeEfficiency;
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
	getImplantedFlux(std::vector<IdType> map)
	{
		std::vector<std::pair<IdType, double>> toReturn;
		// Loop on the map
		for (auto i = 0; i < map.size(); i++) {
			// Look for this value in fluxIndices
			for (auto j = 0; j < fluxIndices.size(); j++) {
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
