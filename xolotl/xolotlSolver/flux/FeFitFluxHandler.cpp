#include "FeFitFluxHandler.h"
#include <iostream>
#include <cmath>

using namespace xolotlSolver;

FeFitFluxHandler::FeFitFluxHandler() {

}

void FeFitFluxHandler::initializeFluxHandler(int numGridpoints, double step) {

	// Set the step size
	stepSize = step;

	// Single He production rate = 0.0029593 * incidentFlux  ( in 1 / nm^3 /sec)
	// where incidentFlux is shown below by Single He flux
	//
	// Single He flux
	// calculated at each depth (x <= 128nm, if x > 128nm, no production)
	double a0 = -0.00073090;
	double a1 = -0.0029330;
	double b1 = 0.0039810;
	double a2 = -0.0041960;
	double b2 = 0.0084180;
	double a3 = -0.0015640;
	double b3 = 0.0099430;
	double a4 = 0.0025910;
	double b4 = 0.0063010;
	double a5 = 0.0038630;
	double b5 = 0.000880;
	double a6 = 0.0022260;
	double b6 = -0.0017580;
	double a7 = 0.00053690;
	double b7 = -0.0013570;
	double a8 = 1.09200000;
	double b8 = -0.00035910;
	double w = 0.013880;

	for (int i = 0; i < numGridpoints; i++) {
		auto x = i * stepSize;
		auto incidentFlux = a0 + a1 * cos(x * w) + b1 * sin(x * w)
				+ a2 * cos(2.0 * x * w) + b2 * sin(2.0 * x * w)
				+ a3 * cos(3.0 * x * w) + b3 * sin(3.0 * x * w)
				+ a4 * cos(4.0 * x * w) + b4 * sin(4.0 * x * w)
				+ a5 * cos(5.0 * x * w) + b5 * sin(5.0 * x * w)
				+ a6 * cos(6.0 * x * w) + b6 * sin(6.0 * x * w)
				+ a7 * cos(7.0 * x * w) + b7 * sin(7.0 * x * w)
				+ a8 * cos(8.0 * x * w) + b8 * sin(8.0 * x * w);
		incidentFluxVec.push_back(incidentFlux);

	}

}
