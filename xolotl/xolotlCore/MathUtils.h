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

#include <array>
#include <limits>
#include <cmath>
#include <numeric>

#include <Kokkos_Array.hpp>

namespace xolotlCore {

//@{
/**
 * Duplicate of the generic std::min and std::max to be used in device kernels
 */
template <typename T>
KOKKOS_INLINE_FUNCTION
const T&
min(const T& a, const T& b)
{
    if (b < a) {
        return b;
    }
    return a;
}

template <typename T>
KOKKOS_INLINE_FUNCTION
const T&
max(const T& a, const T& b)
{
    if (a < b) {
        return b;
    }
    return a;
}
//@}

/**
 * This function tests two doubles to see if they are equal.
 * @param a The first double
 * @param b The second double
 * @return True if the doubles are equal to within machine precision, false
 * otherwise.
 *
 * TODO be careful - this test is effective only if a and b are small
 * (e.g., less than 1).  It is not an effective test if values are large.
 * A relative error test would be more effective.
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
KOKKOS_INLINE_FUNCTION double legendrePolynomial(double x, int degree) {
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

//@{
/**
 * This operation computes the Nth order Legendre polynomials
 *
 * f(x) = c0*P_0(x) + c1*P_1(x) + ... + cN*P_N(x)
 *
 * for a coefficient set {c0,c1,...,cN}.
 *
 * @param x
 *            The x value of the function.
 * @param coeffs
 *            The coefficients array.
 */
template<uint32_t N>
inline double computeNthOrderLegendre(double x,
		const std::array<double, N + 1>& coeffs) {
	int currDegree = 0;
	auto valAtX =
			std::accumulate(coeffs.begin(), coeffs.end(), 0.0,
					[x,&currDegree](double running, double currCoeff) {
						return running + (currCoeff * legendrePolynomial(x, currDegree++));
					});
	return valAtX;
}

template <std::size_t N>
KOKKOS_INLINE_FUNCTION
double
computeNthOrderLegendre(double x, const Kokkos::Array<double, N+1>& coeffs)
{
    int currDegree = 0;
    double value = 0.0;
    for (auto currCoeff : coeffs) {
        value += currCoeff * legendrePolynomial(x, currDegree);
        ++currDegree;
    }
    return value;
}
//@}

/**
 * Computes
 *
 * sum (n - mean) from alpha to beta
 */
KOKKOS_INLINE_FUNCTION
double firstOrderSum(double alpha, double beta, double mean) {
	if (alpha > beta) return 0.0;
	double toReturn = ((beta * (beta + 1.0)) / 2.0)
			- ((alpha * (alpha - 1.0)) / 2.0) - ((beta - alpha + 1.0) * mean);

	return toReturn;
}

/**
 * Computes
 *
 * sum (n - mean)^2 from alpha to beta
 */
KOKKOS_INLINE_FUNCTION
double secondOrderSum(double alpha, double beta, double mean) {
	if (alpha > beta) return 0.0;
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
KOKKOS_INLINE_FUNCTION
double secondOrderOffsetSum(double alpha, double beta, double mean1,
		double mean2, double offset) {
	if (alpha > beta) return 0.0;
	double toReturn = (beta * (beta + 1.0) * (2.0 * beta + 1.0)) / 6.0;

	toReturn -= (alpha * (alpha - 1.0) * (2.0 * alpha - 1.0)) / 6.0;

	toReturn += (offset - mean1 - mean2)
			* ((((beta * (beta + 1.0))) / 2.0)
					- (((alpha * (alpha - 1.0))) / 2.0));

	toReturn += mean1 * (mean2 - offset) * (beta - alpha + 1.0);

	return toReturn;
}

}

#endif /* MATHUTILS_H_ */
