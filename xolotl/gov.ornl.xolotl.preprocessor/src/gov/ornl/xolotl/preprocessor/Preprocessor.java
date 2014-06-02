/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

import java.util.ArrayList;

/**
 * This class generates a valid Xolotl input file with each line representing a
 * cluster with the following variables:
 * 
 * nHe nV nI E_He E_V E_I E_migration D_0
 * 
 * Each element is:
 * 
 * nHe - The number of Helium atoms in the cluster.
 * 
 * nV - The number of vacancies in the cluster.
 * 
 * nI - The number of interstitials in the cluster.
 * 
 * E_He - The binding energy of the cluster with He_1.
 * 
 * E_V - The binding energy of the cluster with V_1.
 * 
 * E_I - The binding energy of the cluster with I_1.
 * 
 * E_migration - The migration energy with which the cluster diffuses.
 * 
 * D_0 - The diffusion factor that sets the scale of the diffusion. (Used to
 * calculate the diffusion coefficient.)
 * 
 * @author Jay Jay Billings
 * 
 */
public class Preprocessor {

	// The maximum size of a Helium cluster in the network.
	private int maxHe = 8;

	// The maximum size of a vacancy cluster in the network.
	private int maxV = 29;

	// The maximum size of an interstitial cluster in the network.
	private int maxI = 2;

	/**
	 * The maximum number of Helium atoms that can be combined with a vacancy
	 * cluster with size equal to the index i in the array. For example, an HeV
	 * size cluster with size 1 would have size = i = 1 and could support a
	 * mixture of up to nine Helium atoms with one vacancy.
	 * 
	 * The special case size = i = 0 is just the maximum size of single Helium.
	 * It is stored in the array for efficiency since it allows us to combine
	 * the loops over helium, vacancies and HeV clusters.
	 */
	private int[] maxHePerV = { maxHe, 9, 14, 18, 20, 27, 30, 35, 40, 45, 50,
			55, 60, 65, 70, 75, 80, 85, 90, 95, 98, 100, 101, 103, 105, 107,
			109, 110, 112, 116 };

	/**
	 * The binding energy engine used to generate binding energies for the
	 * different clusters.
	 */
	private BindingEnergyEngine bindingEnergyEngine = new BindingEnergyEngine();

	/**
	 * This operation generates all He, V and HeV clusters in the network.
	 */
	private ArrayList<Cluster> generateHeV() {

		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Loop over vacancies in the outer loop. Start at zero to account for
		// single He.
		for (int i = 0; i < maxV; ++i) {
			// Add Helium
			for (int j = 0; j < maxHePerV[i]; j++) {
				// Create the cluster
				Cluster tmpCluster = new Cluster();
				tmpCluster.nHe = j + 1;
				tmpCluster.nV = i;
				// Treat everything like a mixed cluster and let the
				// BindingEnergyEngine delegate for single species clusters.
				tmpCluster.E_He = bindingEnergyEngine.getHeVtoHeBindingEnergy(
						j + 1, i);
				tmpCluster.E_V = bindingEnergyEngine.getHeVtoVBindingEnergy(
						j + 1, i);
				// He, V and HeV cannot dissociate into I by normal means.
				tmpCluster.E_I = Double.POSITIVE_INFINITY;
				tmpCluster.E_m = Double.POSITIVE_INFINITY; // FIXME! Need to add
															// migration
															// energies.
				tmpCluster.D_0 = 0.0; // FIXME! Need to add diffusion factors.
				// Add the cluster to the list
				clusterList.add(tmpCluster);
			}
		}

		return clusterList;
	};

	/**
	 * This operation generates all interstitial clusters in the network.
	 */
	private ArrayList<Cluster> generateInterstitials() {

		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Add the interstitials to the list
		for (int i = 0; i < maxI; ++i) {
			// Create the cluster
			Cluster tmpCluster = new Cluster();
			tmpCluster.nI = i + 1;
			tmpCluster.E_I = bindingEnergyEngine.getIBindingEnergy(i + 1);
			// Interstitials cannot dissociate into He and V!
			tmpCluster.E_He = Double.POSITIVE_INFINITY;
			tmpCluster.E_V = Double.POSITIVE_INFINITY;
			tmpCluster.E_m = Double.POSITIVE_INFINITY; // FIXME! Need to add
														// migration energies.
			tmpCluster.D_0 = 0.0; // FIXME! Need to add diffusion factors.
			// Add the cluster to the list
			clusterList.add(tmpCluster);
		}

		return clusterList;
	};

	/**
	 * This operation generates the initial conditions based on the defaults and
	 * the incoming command line arguments.
	 * 
	 * @param args
	 *            The arguments read from the command line
	 * @return The list of clusters created by the preprocessor based on its
	 *         arguments and settings.
	 */
	public ArrayList<Cluster> generate(String[] args) {

		// Create the list of clusters
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Add all of the clusters to the list
		clusterList.addAll(generateInterstitials());
		clusterList.addAll(generateHeV());

		return clusterList;
	}

}
