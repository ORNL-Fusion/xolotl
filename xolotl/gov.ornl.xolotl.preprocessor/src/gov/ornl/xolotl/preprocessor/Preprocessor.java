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
 * All energies are in eV and all diffusion factors are in nm^2/s.
 * 
 * The preprocessor generates 2067 clusters with its default configuration. If
 * only the maximum number of vacancies is changed it generates
 * 
 * maxI + maxHe + 2053 + (maxV - 29) + sum{30<= i <= maxV}(4*i)
 * 
 * clusters.
 * 
 * @author Jay Jay Billings
 * 
 */
public class Preprocessor {

	// The maximum size of a Helium cluster in the network.
	private int maxHe = 8;

	// The maximum size of a mobile He cluster.
	private int maxHeDiffusionSize = 6;

	// The diffusion factors for single species Helium clusters.
	private double[] heDiffusionFactors = { 0.0, 29.0, 32.0, 23.0, 17.0, 5.0,
			10.0 };

	// The migration energies for single species Helium clusters.
	private double[] heMigrationEnergies = { Double.POSITIVE_INFINITY, 0.13,
			0.20, 0.25, 0.20, 0.12, 0.3 };

	// The maximum size of a vacancy cluster in the network.
	private int maxV = 96;

	// The diffusion factor for a single vacancy.
	private double vOneDiffusionFactor = 1800.0;

	// The migration energy for a single vacancy.
	private double vOneMigrationEnergy = 1.30;

	// The maximum size of an interstitial cluster in the network.
	private int maxI = 6;

	// The maximum size of a mobile interstitial cluster.
	private int maxIDiffusionSize = 5;

	// The diffusion coefficients for single species interstitial clusters.
	private double[] iDiffusionFactors = { 0.0, 88.0, 80.0, 39.0, 20.0, 10.0 };

	// The migration energies for single species interstitial clusters.
	private double[] iMigrationEnergies = { Double.POSITIVE_INFINITY, 0.01,
			0.02, 0.03, 0.04, 0.05 };

	/**
	 * The maximum number of Helium atoms that can be combined with a vacancy
	 * cluster with size equal to the index i in the array plus one. For
	 * example, an HeV size cluster with size 1 would have size = i+1 = 1 and i
	 * = 0. It could support a mixture of up to nine Helium atoms with one
	 * vacancy.
	 */
	private int[] maxHePerV = { 9, 14, 18, 20, 27, 30, 35, 40, 45, 50, 55, 60,
			65, 70, 75, 80, 85, 90, 95, 98, 100, 101, 103, 105, 107, 109, 110,
			112, 116 };

	/**
	 * The binding energy engine used to generate binding energies for the
	 * different clusters.
	 */
	private BindingEnergyEngine bindingEnergyEngine = new BindingEnergyEngine();

	/**
	 * This operation generates all Helium clusters in the network.
	 * 
	 * @return A list of clusters configured to satisfy the bounds and composed
	 *         solely of Helium.
	 */
	private ArrayList<Cluster> generateHe() {

		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Create the He clusters
		for (int i = 0; i < maxHe; i++) {
			// Create the cluster
			Cluster tmpCluster = new Cluster();
			tmpCluster.nHe = i + 1;
			tmpCluster.E_He = bindingEnergyEngine.getHeBindingEnergy(i + 1);
			// Add the cluster to the list
			clusterList.add(tmpCluster);
		}

		// Configure the diffusion parameters.
		for (int i = 0; i < maxHeDiffusionSize; i++) {
			Cluster tmpCluster = clusterList.get(i);
			tmpCluster.D_0 = heDiffusionFactors[i + 1];
			tmpCluster.E_m = heMigrationEnergies[i + 1];
		}

		return clusterList;
	}

	/**
	 * This operation creates an HeV cluster with the specified size. It
	 * configures the binding energies on its own.
	 * 
	 * @param heSize
	 *            The number of Helium atoms in the cluster
	 * @param vSize
	 *            The number of vacancies in the cluster
	 * @return The cluster.
	 */
	private Cluster makeHeVCluster(int heSize, int vSize) {
		// Create the cluster
		Cluster cluster = new Cluster();
		cluster.nHe = heSize;
		cluster.nV = vSize;
		// Treat everything like a mixed cluster and let the
		// BindingEnergyEngine delegate for single species clusters.
		cluster.E_He = bindingEnergyEngine.getHeVtoHeBindingEnergy(heSize,
				vSize);
		cluster.E_V = bindingEnergyEngine.getHeVtoVBindingEnergy(heSize, vSize);

		return cluster;
	}

	/**
	 * This operation generates all V and HeV clusters in the network.
	 * 
	 * @return A list of clusters configured to satisfy the bounds and composed
	 *         solely of V and HeV clusters.
	 */
	private ArrayList<Cluster> generateHeV() {

		// Local Declarations
		ArrayList<Cluster> clusterList = new ArrayList<Cluster>();

		// Loop over vacancies in the outer loop. Start at zero to account for
		// single He. This creates V and HeV up to the maximum size in the
		// maxHePerV array.
		for (int i = 1; i <= maxV && i <= maxHePerV.length; ++i) {
			// Add Helium
			for (int j = 0; j <= maxHePerV[i - 1]; j++) {
				// Add the cluster to the list
				clusterList.add(makeHeVCluster(j, i));
			}
		}

		// Create V and HeV up to the maximum length with a constant nHe/nV = 4.
		for (int i = maxHePerV.length + 1; i <= maxV; i++) {
			// Add Helium
			for (int j = 0; j <= i * 4; j++) {
				// Add the cluster to the list
				clusterList.add(makeHeVCluster(j, i));
			}
		}

		// Set V_1 diffusion parameters. V_1 is the first in the list, so it is
		// straightforward to set it.
		clusterList.get(0).D_0 = vOneDiffusionFactor;
		clusterList.get(0).E_m = vOneMigrationEnergy;

		return clusterList;
	};

	/**
	 * This operation generates all interstitials clusters in the network.
	 * 
	 * @return A list of clusters configured to satisfy the bounds and composed
	 *         solely of interstitials.
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
			// Add the cluster to the list
			clusterList.add(tmpCluster);
		}

		// Configure the diffusion parameters.
		for (int i = 0; i < maxIDiffusionSize; i++) {
			Cluster tmpCluster = clusterList.get(i);
			tmpCluster.D_0 = iDiffusionFactors[i + 1];
			tmpCluster.E_m = iMigrationEnergies[i + 1];
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
		clusterList.addAll(generateHe());
		clusterList.addAll(generateHeV());

		return clusterList;
	}

}
