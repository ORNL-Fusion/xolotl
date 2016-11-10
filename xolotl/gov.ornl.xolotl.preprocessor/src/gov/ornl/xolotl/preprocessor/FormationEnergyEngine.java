/**
 * 
 */
package gov.ornl.xolotl.preprocessor;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.io.FileInputStream;
import java.io.IOException;

/**
 * This class computes the formation energies for a cluster with a given
 * composition. All formation energies are in electron volts (eV).
 * 
 * This class looks in the user directory during construction and will override
 * the fit coefficients if it finds a csv file with one set of coefficients per
 * row and named fit.csv. There should be eight rows with six coefficients per
 * row.
 * 
 * @author Jay Jay Billings
 * 
 */
public class FormationEnergyEngine {

	/**
	 * The set of Helium formation energies up to He_8 indexed by size. That is
	 * E_(f,He_1) = heFormationEnergies[1]. The value at index zero is just
	 * padding to make the indexing easy.
	 */
	private static double[] heFormationEnergies = { Double.POSITIVE_INFINITY, 6.15, 11.44, 16.35, 21.0, 26.1, 30.24,
			34.93, 38.80 };

	/**
	 * The set of xenon formation energies up to Xe_29 indexed by size. That is
	 * E_(f,Xe_1) = xeFormationEnergies[1]. The value at index zero is just
	 * padding to make the indexing easy.
	 */
	private static double[] xeFormationEnergies = { Double.POSITIVE_INFINITY, 7.0, 12.15, 17.15, 21.90, 26.50, 31.05,
			35.30, 39.45, 43.00, 46.90, 50.65, 53.90, 56.90, 59.80, 62.55, 65.05, 67.45, 69.45, 71.20, 72.75, 74.15,
			75.35, 76.40, 77.25, 77.95, 78.45, 78.80, 78.95, 79.0 };

	/**
	 * The set of vacancy formation energies up to V_2 indexed by size. That is
	 * E_(f,V_1) = vFormationEnergies[1]. The value at index zero is just
	 * padding to make the indexing easy.
	 */
	private static double[] vFormationEnergies = { Double.POSITIVE_INFINITY, 3.6, 7.25 };

	/**
	 * The set of interstitial formation energies up to I_6 indexed by size.
	 * That is E_(f,V_1) = iFormationEnergies[1]. The value at index zero is
	 * just padding to make the indexing easy.
	 */
	private static double[] iFormationEnergies = { Double.POSITIVE_INFINITY, 10.0, 18.5, 27.0, 35.0, 42.5, 48.0 };

	// Coefficients for the Legendre polynomial fit
	// Low means V <= 27
	// Coefficients for c_0 in the 2D E_f,HeV fit
	private double[] c0CoefficientsLow = { 253.35, 435.36, 336.50, 198.92, 95.154, 21.544 };
	// Coefficients for c_1 in the 2D E_f,HeV fit
	private double[] c1CoefficientsLow = { 493.29, 1061.3, 1023.9, 662.92, 294.24, 66.962 };
	// Coefficients for c_2 in the 2D E_f,HeV fit
	private double[] c2CoefficientsLow = { 410.40, 994.89, 1044.6, 689.41, 286.52, 60.712 };
	// Coefficients for c_3 in the 2D E_f,HeV fit
	private double[] c3CoefficientsLow = { 152.99, 353.16, 356.10, 225.75, 87.077, 15.640 };
	// High means V > 27
	// Coefficients for c_0 in the 2D E_f,HeV fit
	private double[] c0CoefficientsHigh = { -847.90, -3346.9, -4510.3, -3094.7, -971.18, -83.770 };
	// Coefficients for c_1 in the 2D E_f,HeV fit
	private double[] c1CoefficientsHigh = { -1589.3, -4894.6, -6001.8, -4057.5, -1376.4, -161.91 };
	// Coefficients for c_2 in the 2D E_f,HeV fit
	private double[] c2CoefficientsHigh = { 834.91, 1981.8, 1885.7, 1027.1, 296.69, 29.902 };
	// Coefficients for c_3 in the 2D E_f,HeV fit
	private double[] c3CoefficientsHigh = { 1547.2, 3532.3, 3383.6, 1969.2, 695.17, 119.23 };

	/**
	 * The formation energies for He_xV_1. The entry at i = 0 is for a single
	 * vacancy (He_0V_1) and is there as a buffer. Like the formation energies,
	 * i = heSize.
	 */
	private double[] heV1FormationEnergies = { vFormationEnergies[1], 5.14166, 8.20919, 11.5304, 14.8829, 18.6971,
			22.2847, 26.3631, 30.1049, 34.0081, 38.2069, 42.4217, 46.7378, 51.1551, 55.6738 };

	/**
	 * The formation energies for He_xV_2. The entry at i = 0 is for a
	 * di-vacancy (He_0V_2) and is there as a buffer. Like the formation
	 * energies, i = heSize.
	 */
	private double[] heV2FormationEnergies = { vFormationEnergies[2], 7.10098, 8.39913, 9.41133, 11.8748, 14.8296,
			17.7259, 20.7747, 23.7993, 26.7984, 30.0626, 33.0385, 36.5173, 39.9406, 43.48, 46.8537, 50.4484, 54.0879,
			57.7939 };

	/**
	 * The constructor. If the file fit.csv exists in user.dir, this operation
	 * will try to load the coefficients from it. It closes the file when it is
	 * finished.
	 */
	public FormationEnergyEngine() {
		// Get the fit file
		File fitFile = new File("fit.csv");

		// Load the file if it exists
		if (fitFile.exists()) {
			try {
				// Create the readers
				InputStreamReader fitStreamReader = new InputStreamReader(new FileInputStream(fitFile), "UTF-8");
				BufferedReader fitReader = new BufferedReader(fitStreamReader);

				// Read from the Legendre fit file
				String[] c0Strings = fitReader.readLine().split(",");
				String[] c1Strings = fitReader.readLine().split(",");
				String[] c2Strings = fitReader.readLine().split(",");
				String[] c3Strings = fitReader.readLine().split(",");
				String[] c4Strings = fitReader.readLine().split(",");
				String[] c5Strings = fitReader.readLine().split(",");
				String[] c6Strings = fitReader.readLine().split(",");
				String[] c7Strings = fitReader.readLine().split(",");

				// Convert the c0 coefficients
				c0CoefficientsLow[0] = Double.valueOf(c0Strings[0]);
				c0CoefficientsLow[1] = Double.valueOf(c0Strings[1]);
				c0CoefficientsLow[2] = Double.valueOf(c0Strings[2]);
				c0CoefficientsLow[3] = Double.valueOf(c0Strings[3]);
				c0CoefficientsLow[4] = Double.valueOf(c0Strings[4]);
				c0CoefficientsLow[5] = Double.valueOf(c0Strings[5]);
				// Convert the c1 coefficients
				c1CoefficientsLow[0] = Double.valueOf(c1Strings[0]);
				c1CoefficientsLow[1] = Double.valueOf(c1Strings[1]);
				c1CoefficientsLow[2] = Double.valueOf(c1Strings[2]);
				c1CoefficientsLow[3] = Double.valueOf(c1Strings[3]);
				c1CoefficientsLow[4] = Double.valueOf(c1Strings[4]);
				c1CoefficientsLow[5] = Double.valueOf(c1Strings[5]);
				// Convert the c2 coefficients
				c2CoefficientsLow[0] = Double.valueOf(c2Strings[0]);
				c2CoefficientsLow[1] = Double.valueOf(c2Strings[1]);
				c2CoefficientsLow[2] = Double.valueOf(c2Strings[2]);
				c2CoefficientsLow[3] = Double.valueOf(c2Strings[3]);
				c2CoefficientsLow[4] = Double.valueOf(c2Strings[4]);
				c2CoefficientsLow[5] = Double.valueOf(c2Strings[5]);
				// Convert the c3 coefficients
				c3CoefficientsLow[0] = Double.valueOf(c3Strings[0]);
				c3CoefficientsLow[1] = Double.valueOf(c3Strings[1]);
				c3CoefficientsLow[2] = Double.valueOf(c3Strings[2]);
				c3CoefficientsLow[3] = Double.valueOf(c3Strings[3]);
				c3CoefficientsLow[4] = Double.valueOf(c3Strings[4]);
				c3CoefficientsLow[5] = Double.valueOf(c3Strings[5]);
				// Convert the c0 coefficients
				c0CoefficientsHigh[0] = Double.valueOf(c4Strings[0]);
				c0CoefficientsHigh[1] = Double.valueOf(c4Strings[1]);
				c0CoefficientsHigh[2] = Double.valueOf(c4Strings[2]);
				c0CoefficientsHigh[3] = Double.valueOf(c4Strings[3]);
				c0CoefficientsHigh[4] = Double.valueOf(c4Strings[4]);
				c0CoefficientsHigh[5] = Double.valueOf(c4Strings[5]);
				// Convert the c1 coefficients
				c1CoefficientsHigh[0] = Double.valueOf(c5Strings[0]);
				c1CoefficientsHigh[1] = Double.valueOf(c5Strings[1]);
				c1CoefficientsHigh[2] = Double.valueOf(c5Strings[2]);
				c1CoefficientsHigh[3] = Double.valueOf(c5Strings[3]);
				c1CoefficientsHigh[4] = Double.valueOf(c5Strings[4]);
				c1CoefficientsHigh[5] = Double.valueOf(c5Strings[5]);
				// Convert the c2 coefficients
				c2CoefficientsHigh[0] = Double.valueOf(c6Strings[0]);
				c2CoefficientsHigh[1] = Double.valueOf(c6Strings[1]);
				c2CoefficientsHigh[2] = Double.valueOf(c6Strings[2]);
				c2CoefficientsHigh[3] = Double.valueOf(c6Strings[3]);
				c2CoefficientsHigh[4] = Double.valueOf(c6Strings[4]);
				c2CoefficientsHigh[5] = Double.valueOf(c6Strings[5]);
				// Convert the c3 coefficients
				c3CoefficientsHigh[0] = Double.valueOf(c7Strings[0]);
				c3CoefficientsHigh[1] = Double.valueOf(c7Strings[1]);
				c3CoefficientsHigh[2] = Double.valueOf(c7Strings[2]);
				c3CoefficientsHigh[3] = Double.valueOf(c7Strings[3]);
				c3CoefficientsHigh[4] = Double.valueOf(c7Strings[4]);
				c3CoefficientsHigh[5] = Double.valueOf(c7Strings[5]);

				// Close the readers
				fitReader.close();
				fitStreamReader.close();
			} catch (FileNotFoundException e) {
				// Complain
				System.err.println("FormationEnergyEngine Message:"
						+ " There is something wrong with your fit.csv file! " + " It does not exist. Aborting.");
				e.printStackTrace();
			} catch (IOException e) {
				// Complain
				System.err.println("FormationEnergyEngine Message:"
						+ " There is something wrong with your fit.csv file! " + " It cannot be read. Aborting.");
				e.printStackTrace();
			}
		}

		return;
	};

	/**
	 * This operation computes and returns the formation energy of a helium
	 * cluster of the specified size.
	 * 
	 * @param size
	 *            The size of the Helium cluster.
	 * @return The formation energy.
	 */
	public double getHeFormationEnergy(int size) {
		// Initialize the formation energy to infinity
		double energy = Double.POSITIVE_INFINITY;

		if (size < 9 && size > 0)
			energy = heFormationEnergies[size];

		return energy;
	}

	/**
	 * This operation computes and returns the formation energy of a xenon
	 * cluster of the specified size.
	 * 
	 * @param size
	 *            The size of the xenon cluster.
	 * @return The formation energy.
	 */
	public double getXeFormationEnergy(int size) {
		// Initialize the formation energy to infinity
		double energy = Double.POSITIVE_INFINITY;

		if (size < 30 && size > 0)
			energy = xeFormationEnergies[size];
		else
			energy = 79.0;

		return energy;
	}

	/**
	 * This operation computes and returns the formation energy of a vacancy
	 * cluster of the specified size.
	 * 
	 * @param size
	 *            The size of the vacancy cluster.
	 * @return The formation energy.
	 */
	public double getVFormationEnergy(int size) {
		// Initialize the formation energy to infinity
		double energy = Double.POSITIVE_INFINITY;

		if (size < 3 && size > 0)
			energy = vFormationEnergies[size];
		else if (size >= 3) {
			// Get the energy for He_xV_y
			energy = getHeVFormationEnergy(0, size);
		}

		return energy;
	}

	/**
	 * This operation returns the formation energy for a cluster composed of
	 * only interstitial defects.
	 * 
	 * The values are taken from the iFormationEnergies array for a size smaller
	 * or equal to 5, and computed with the formula E_f = 48 + 6 * (size - 6) eV
	 * starting at size = 6.
	 * 
	 * @param size
	 *            the number of interstitials in the cluster
	 * @return the formation energy
	 */
	public double getIFormationEnergy(int size) {
		// Initialize the energy
		double energy = Double.POSITIVE_INFINITY;

		// The formation energies from the table are valid for 0 < size < 6.
		if (size > 0 && size < 6) {
			energy = iFormationEnergies[size];
		} else if (size > 5)
			energy = 48.0 + 6.0 * (size - 6.0);

		return energy;
	}

	/**
	 * This operation computes the Legendre polynomial, P_n (x), of degree n.
	 * 
	 * The Legendre polynomials of degree 0 and 1 are P_0(x)=1.0 and P_1(x)=x,
	 * respectively. With these conditions, the Legendre polynomials satisfy the
	 * following recurrence relation: (n+1)*P_(n+1)(x) = (2n+1)*x*P_n(x) -
	 * n*P_(n-1)(x)
	 * 
	 * @param x
	 *            The x value of the function
	 * @param order
	 *            The order of the polynomial
	 */
	private double legendrePolynomial(double x, int degree) {
		// For degree 0 the return value is 1.0
		if (degree == 0)
			return 1.0;
		// For degree1 the return value i x
		if (degree == 1)
			return x;

		// Initialize the polynomials orders for the loop
		double Pn2 = 1.0, Pn1 = x, Pn = 0.0;
		// Loop on the wanted degree
		for (int n = 1; n < degree; n++) {
			// Compute the polynomial at the current order
			Pn = (((2.0 * (double) n + 1.0) * x * Pn1) - ((double) n * Pn2)) / ((double) n + 1.0);
			// Update the polynomials orders
			Pn2 = Pn1;
			Pn1 = Pn;
		}

		return Pn;
	}

	/**
	 * This operation computes the 3rd order Legendre polynomials
	 * 
	 * f(x) = c0*P_0(x) + c1*P_1(x) + c2*P_2(x) + c3*P_3(x) = c0 + c1 * x +
	 * c2*P_2(x) + c3*P_3(x)
	 * 
	 * for a coefficient set {c0,c1,c2,c3}.
	 * 
	 * @param x
	 *            The x value of the function.
	 * @param coeffs
	 *            The coefficients array.
	 */
	private double compute3rdOrderLegendre(double x, double[] coeffs) {
		// Initialize the value
		double value = 0.0;

		// Compute the value
		for (int i = 0; i < 4; i++) {
			value = value + coeffs[i] * legendrePolynomial(x, i);
		}

		return value;
	}

	/**
	 * This operation computes the 5th order Legendre polynomials
	 * 
	 * f(x) = c0*P_0(x) + c1*P_1(x) + c2*P_2(x) + c3*P_3(x) + c4*P_4(x) +
	 * c5*P_5(x) = c0 + c1 * x + c2*P_2(x) + c3*P_3(x) + c4*P_4(x) + c5*P_5(x)
	 * 
	 * for a coefficient set {c0,c1,c2,c3,c4,c5}.
	 * 
	 * @param x
	 *            The x value of the function.
	 * @param coeffs
	 *            The coefficients array.
	 */
	private double compute5thOrderLegendre(double x, double[] coeffs) {
		// Initialize the value
		double value = 0.0;

		// Compute the value
		for (int i = 0; i < 6; i++) {
			value = value + coeffs[i] * legendrePolynomial(x, i);
		}

		return value;
	}

	/**
	 * This operation returns the formation energy for an HeV cluster.
	 * 
	 * The formation energies are determined using a 2D Legendre polynomial fit
	 * f(x/y,y) where y = vSize and x/y = heSize/vSize.
	 * 
	 * @param heSize
	 *            The number of helium atoms in the cluster
	 * @param vSize
	 *            The number of vacancies in the cluster
	 * @return The formation energy
	 */
	public double getHeVFormationEnergy(int heSize, int vSize) {

		double energy = Double.NEGATIVE_INFINITY;
		// The following coefficients are computed using the above and are used
		// to evaluate the full function f(x,y).
		double[] coefficients = { 0.0, 0.0, 0.0, 0.0 };

		// Check to see if the vacancy size is large enough that the energy can
		// be computed from the fit or if it is so small that the exact values
		// must be used instead.
		if (vSize > 2) {
			// Get the He/V ratio
			double x = 2.0 * (((double) heSize / (double) vSize) / 9.0) - 1.0;
			// Initialize the vacancy number
			double y = 0.0;

			// We have 2 fits, one for low V and one for high V
			if (vSize <= 27) {
				// Compute the vacancy number
				y = 2.0 * (((double) vSize - 1.0) / 26.0) - 1.0;
				// Get the coefficients
				coefficients[0] = compute5thOrderLegendre(x, c0CoefficientsLow);
				coefficients[1] = compute5thOrderLegendre(x, c1CoefficientsLow);
				coefficients[2] = compute5thOrderLegendre(x, c2CoefficientsLow);
				coefficients[3] = compute5thOrderLegendre(x, c3CoefficientsLow);
			} else {
				// Compute the vacancy number
				y = 2.0 * (((double) vSize - 1.0) / 451.0) - 1.0;
				// Get the coefficients
				coefficients[0] = compute5thOrderLegendre(x, c0CoefficientsHigh);
				coefficients[1] = compute5thOrderLegendre(x, c1CoefficientsHigh);
				coefficients[2] = compute5thOrderLegendre(x, c2CoefficientsHigh);
				coefficients[3] = compute5thOrderLegendre(x, c3CoefficientsHigh);
			}
			// Get the energy
			energy = compute3rdOrderLegendre(y, coefficients);

		} else if ((vSize == 1 && heSize < heV1FormationEnergies.length)
				|| (vSize == 2 && heSize < heV2FormationEnergies.length)) {
			// Get the exact energy
			energy = (vSize == 1) ? heV1FormationEnergies[heSize] : heV2FormationEnergies[heSize];
		}

		return energy;
	}

}
