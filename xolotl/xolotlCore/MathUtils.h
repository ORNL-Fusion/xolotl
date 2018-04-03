/*
 * MathUtils.h
 *
 * Various math utilities.
 *
 *  Created on: May 20, 2014
 *      Author: Jay Jay Billings
 */

#ifndef MATHUTILS_H_
#define MATHUTILS_H_

#include <limits>
#include <cmath>

namespace xolotlCore {

/**
 * This function tests two doubles to see if they are equal.
 * @param a The first double
 * @param b The second double
 * @return True if the doubles are equal to within machine precision, false
 * otherwise.
 */
inline bool equal(double a, double b) {
	return std::fabs(b - a) < std::numeric_limits<double>::epsilon();
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
inline double legendrePolynomial(double x, int degree) {
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
		Pn = (((2.0 * (double) n + 1.0) * x * Pn1) - ((double) n * Pn2))
				/ ((double) n + 1.0);
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
inline double compute3rdOrderLegendre(double x, std::vector<double> coeffs) {
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
inline double compute5thOrderLegendre(double x, std::vector<double> coeffs) {
	// Initialize the value
	double value = 0.0;

	// Compute the value
	for (int i = 0; i < 6; i++) {
		value = value + coeffs[i] * legendrePolynomial(x, i);
	}

	return value;
}

/**
 * Computes
 *
 * sum (n - mean) from alpha to beta
 */
inline double firstOrderSum(double alpha, double beta, double mean) {
	double toReturn = ((beta * (beta + 1.0)) / 2.0)
			- ((alpha * (alpha - 1.0)) / 2.0) - ((beta - alpha + 1.0) * mean);

	return toReturn;
}

/**
 * Computes
 *
 * sum (n - mean)^2 from alpha to beta
 */
inline double secondOrderSum(double alpha, double beta, double mean) {
	double toReturn = (beta * (beta + 1.0) * ((2.0 * beta) + 1.0)) / 6.0;

	toReturn -= (alpha * (alpha - 1.0) * ((2.0 * alpha) - 1.0)) / 6.0;

	toReturn -= 2.0 * mean
			* ((((beta * (beta + 1.0))) / 2.0)
					- (((alpha * (alpha - 1.0))) / 2.0));

	toReturn += mean * mean * (beta - alpha + 1.0);

	return toReturn;
}

/**
 * Computes
 *
 * sum (n - mean1) * (n + offset - mean2) from alpha to beta
 */
inline double secondOrderOffsetSum(double alpha, double beta, double mean1,
		double mean2, double offset) {
	double toReturn = (beta * (beta + 1.0) * (2.0 * beta + 1.0)) / 6.0;

	toReturn -= (alpha * (alpha - 1.0) * (2.0 * alpha - 1.0)) / 6.0;

	toReturn += (offset - mean1 - mean2)
			* ((((beta * (beta + 1.0))) / 2.0)
					- (((alpha * (alpha - 1.0))) / 2.0));

	toReturn += mean1 * (mean2 - offset) * (beta - alpha + 1.0);

	return toReturn;
}

// Concise names for multi-dimensional arrays.
// For dimensions > 1, this gives C-style row-major ordering.
// (I.e., last dimension varies fastest when indexing.)
template<class T, uint32_t N>
using Array1D = std::array<T, N>;

template<class T, uint32_t N0, uint32_t N1>
using Array2D = Array1D<Array1D<T, N1>, N0>;

template<class T, uint32_t N0, uint32_t N1, uint32_t N2>
using Array3D = Array1D<Array2D<T, N1, N2> , N0>;

}

#endif /* MATHUTILS_H_ */
