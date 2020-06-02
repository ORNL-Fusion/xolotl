#ifndef OPTIONS_H
#define OPTIONS_H

#include <xolotl/options/IOptions.h>

namespace xolotl {
namespace options {

/**
 * Options realizes the IOptions interface.
 * All private members will be accessed through getters and setters.
 */
class Options: public IOptions {

protected:

	/**
	 * The flag that says if Xolotl should run.
	 */
	bool shouldRunFlag;

	/**
	 * The value of the exit code. Should be 0 if everything went well.
	 */
	int exitCode;

	/**
	 * The name of the file where the network is stored.
	 */
	std::string networkFilename;

	/**
	 * The options that will be given to PETSc.
	 */
	std::string petscArg;

	/**
	 * Use the constant temperature set of handlers?
	 */
	bool constTempFlag;

	/**
	 * Value for the constant temperature.
	 */
	double constTemperature;

	/**
	 * Use the temperature profile set of handlers?
	 */
	bool tempProfileFlag;

	/**
	 * Name of the input temperature profile file.
	 */
	std::string tempProfileFilename;

	/**
	 * Use the heat equation set of handlers?
	 */
	bool heatFlag;

	/**
	 * Value for the bulk temperature.
	 */
	double bulkTemperature;

	/**
	 * Use the flux amplitude option?
	 */
	bool fluxFlag;

	/**
	 * Value for the  flux.
	 */
	double fluxAmplitude;

	/**
	 * Use a time profile for the flux?
	 */
	bool fluxProfileFlag;

	/**
	 * Name of the input time profile file for the flux.
	 */
	std::string fluxProfileFilename;

	/**
	 * Which type of performance infrastructure should we use?
	 */
	perf::IHandlerRegistry::RegistryType perfRegistryType;

	/**
	 * Use the "standard" set of handlers for the visualization infrastructure?
	 */
	bool vizStandardHandlersFlag;

	/**
	 * Name of the material.
	 */
	std::string materialName;

	/**
	 * Value of the initial vacancy concentration.
	 */
	double initialVConcentration;

	/**
	 * Value of the electronic stopping power.
	 */
	double zeta;

	/**
	 * Value of the portion of the void on the grid at the start of the simulation.
	 */
	double voidPortion;

	/**
	 * Number of dimensions for the simulation.
	 */
	int dimensionNumber;

	/**
	 * Use a regular grid on the x direction?
	 */
	bool useRegularGridFlag;

	/**
	 * Use a Chebyshev grid on the x direction?
	 */
	bool useChebyshevGridFlag;

	/**
	 * Read in the grid on the x direction?
	 */
	bool readInGridFlag;

	/**
	 * The name of the file where the grid is stored.
	 */
	std::string gridFilename;

	/**
	 * The map of physical processes to use in the simulation.
	 */
	std::map<std::string, bool> processMap;

	/**
	 * String of the list of wanted GB.
	 */
	std::string gbList;

	/**
	 * Minimum size for the grouping.
	 */
	int groupingMin;

	/**
	 * Width for the grouping in the first direction.
	 */
	int groupingWidthA;

	/**
	 * Width for the grouping in the second direction.
	 */
	int groupingWidthB;

	/**
	 * Value of the sputtering yield.
	 */
	double sputteringYield;

	/**
	 * Use a HDF5 file?
	 */
	bool useHDF5Flag;

	/**
	 * Use the phase cut for the network?
	 */
	bool usePhaseCutFlag;

	/**
	 * Maximum number of He or Xe
	 */
	int maxImpurity;

	/**
	 * Maximum number of D
	 */
	int maxD;

	/**
	 * Maximum number of T
	 */
	int maxT;

	/**
	 * Maximum number of V
	 */
	int maxV;

	/**
	 * Maximum number of I
	 */
	int maxI;

	/**
	 * Number of grid point in the depth direction
	 */
	int nX;

	/**
	 * Step size in the depth direction
	 */
	double xStepSize;

	/**
	 * Number of grid point in the Y direction
	 */
	int nY;

	/**
	 * Step size in the Y direction
	 */
	double yStepSize;

	/**
	 * Number of grid point in the Z direction
	 */
	int nZ;

	/**
	 * Step size in the Z direction
	 */
	double zStepSize;

	/**
	 * The boundary condition on the given side of the grid.
	 */
	int leftBoundary, rightBoundary, bottomBoundary, topBoundary, frontBoundary,
			backBoundary;

	/**
	 * Depth for the bubble bursting in nm.
	 */
	double burstingDepth;

	/**
	 * An explicitly-given value to use to seed the random number generator.
	 * Only used if rngUseSeed is true.
	 */
	unsigned int rngSeed;

	/**
	 * Whether to use the rngSeed value to seed the random number generator.
	 */
	bool rngUseSeed;

	/**
	 * Whether to print the value used to seed the random number generator.
	 */
	bool rngPrintSeed;

	/**
	 * Re-solution minimum size
	 */
	int resoMinSize;

	/**
	 * Average radius computation minimum size
	 */
    util::Array<int, 4> radiusMinSizes;

	/**
	 * Density of atom in a bubble in nm-3.
	 */
	double density;

	/**
	 * Length of time of the pulse in s.
	 */
	double pulseTime;

	/**
	 * Proportion of the pulse that is on.
	 */
	double pulseProportion;

	/**
	 * Length of the lattice side in nm.
	 */
	double latticeParameter;

	/**
	 * Radius of the main impurity (He, Xe) in nm.
	 */
	double impurityRadius;

	/**
	 * Reflect the fact that interstitial clusters have a larger surrounding strain field.
	 */
	double biasFactor;

	/**
	 * Factor between the He and the H radius.
	 */
	double hydrogenFactor;

	/**
	 * Xenon diffusion coefficient in nm2 s-1
	 */
	double xenonDiffusivity;

	/**
	 * Fission yield, how many xenon atoms are created per fission
	 */
	double fissionYield;

	/**
	 * Migration energy above which the diffusion will be ignored
	 */
	double migrationThreshold;

public:

	/**
	 * The constructor.
	 */
	Options();

	/**
	 * The destructor.
	 */
	~Options();

	/**
	 * Read the parameters from the given file to set the different
	 * xolotl options.
	 * \see IOptions.h
	 */
	void readParams(int argc, char *argv[]) override;

	/**
	 * Should the program run after parsing the parameter file?
	 * \see IOptions.h
	 */
	bool shouldRun() const override {
		return shouldRunFlag;
	}

	/**
	 * Set the shouldRunFlag.
	 * \see IOptions.h
	 */
	void setShouldRunFlag(bool flag) override {
		shouldRunFlag = flag;
	}

	/**
	 * If program shouldn't run, what should its exit code be?
	 * \see IOptions.h
	 */
	int getExitCode() const override {
		return exitCode;
	}

	/**
	 * Set the value for the exit code.
	 * \see IOptions.h
	 */
	void setExitCode(int code) override {
		exitCode = code;
	}

	/**
	 * Get the name of the network file.
	 * \see IOptions.h
	 */
	std::string getNetworkFilename() const override {
		return networkFilename;
	}

	/**
	 * Get the Arguments for PETSc.
	 * \see IOptions.h
	 */
	std::string getPetscArg() const override {
		return petscArg;
	}

	/**
	 * Should we use const temperature handlers?
	 * \see IOptions.h
	 */
	bool useConstTemperatureHandlers() const override {
		return constTempFlag;
	}

	/**
	 * Obtain the value of the constant temperature to be used.
	 * \see IOptions.h
	 */
	double getConstTemperature() const override {
		return constTemperature;
	}

	/**
	 * Should we use temperature profile handlers?
	 * \see IOptions.h
	 */
	bool useTemperatureProfileHandlers() const override {
		return tempProfileFlag;
	}

	/**
	 * Obtain the name of the file containing the temperature profile data.
	 * \see IOptions.h
	 */
	std::string getTempProfileFilename() const override {
		return tempProfileFilename;
	}

	/**
	 * Should we use heat equation handlers?
	 * \see IOptions.h
	 */
	bool useHeatEquationHandlers() const override {
		return heatFlag;
	}

	/**
	 * Obtain the value of the temperature to be used in the bulk.
	 * \see IOptions.h
	 */
	double getBulkTemperature() const override {
		return bulkTemperature;
	}

	/**
	 * Should we use the flux option?
	 * \see IOptions.h
	 */
	bool useFluxAmplitude() const override {
		return fluxFlag;
	}

	/**
	 * Obtain the value of the flux intensity to be used.
	 * \see IOptions.h
	 */
	double getFluxAmplitude() const override {
		return fluxAmplitude;
	}

	/**
	 * Should we use a time profile for the flux?
	 * \see IOptions.h
	 */
	bool useFluxTimeProfile() const override {
		return fluxProfileFlag;
	}

	/**
	 * Obtain the name of the file containing the time profile data for the
	 * flux.
	 * \see IOptions.h
	 */
	std::string getFluxProfileName() const override {
		return fluxProfileFilename;
	}

	/**
	 * Which type of performance handlers should we use?
	 * \see IOptions.h
	 */
	perf::IHandlerRegistry::RegistryType getPerfHandlerType(void) const
			override {
		return perfRegistryType;
	}

	/**
	 * Should we use the "standard" set of handlers for the visualization?
	 * If false, use dummy (stub) handlers.
	 * \see IOptions.h
	 */
	bool useVizStandardHandlers() const override {
		return vizStandardHandlersFlag;
	}

	/**
	 * Obtain the name of the material to be used for simulation.
	 * \see IOptions.h
	 */
	std::string getMaterial() const override {
		return materialName;
	}

	/**
	 * Obtain the value of the concentration for the vacancies.
	 * \see IOptions.h
	 */
	double getInitialVConcentration() const override {
		return initialVConcentration;
	}

	/**
	 * Obtain the value of the electronic stopping power.
	 * \see IOptions.h
	 */
	double getZeta() const override {
		return zeta;
	}

	/**
	 * Obtain the number of dimensions for the simulation.
	 * \see IOptions.h
	 */
	int getDimensionNumber() const override {
		return dimensionNumber;
	}

	/**
	 * Obtain the value of the void portion for the simulation.
	 * \see IOptions.h
	 */
	double getVoidPortion() const override {
		return voidPortion;
	}

	/**
	 * Should we use a regular grid on the x direction?
	 * \see IOptions.h
	 */
	bool useRegularXGrid() const override {
		return useRegularGridFlag;
	}

	/**
	 * Should we use a Chebyshev grid on the x direction?
	 * \see IOptions.h
	 */
	bool useChebyshevGrid() const override {
		return useChebyshevGridFlag;
	}

	/**
	 * Should we read in the grid on the x direction?
	 * \see IOptions.h
	 */
	bool useReadInGrid() const override {
		return readInGridFlag;
	}

	/**
	 * Get the name of the grid file.
	 * \see IOptions.h
	 */
	std::string getGridFilename() const override {
		return gridFilename;
	}

	/**
	 * Obtain the physical process map.
	 *
	 * @return The map
	 */
	std::map<std::string, bool> getProcesses() const override {
		return processMap;
	}

	/**
	 * Obtain the string listing the wanted GB.
	 * \see IOptions.h
	 */
	std::string getGbString() const override {
		return gbList;
	}

	/**
	 * Obtain the minimum size for the grouping.
	 * \see IOptions.h
	 */
	int getGroupingMin() const override {
		return groupingMin;
	}

	/**
	 * Obtain the first width for the grouping.
	 * \see IOptions.h
	 */
	int getGroupingWidthA() const override {
		return groupingWidthA;
	}

	/**
	 * Obtain the second width for the grouping.
	 * \see IOptions.h
	 */
	int getGroupingWidthB() const override {
		return groupingWidthB;
	}

	/**
	 * Obtain the value of the intensity of the sputtering yield to be used.
	 * \see IOptions.h
	 */
	double getSputteringYield() const override {
		return sputteringYield;
	}

	/**
	 * To know if we should use the HDF5 file.
	 * \see IOptions.h
	 */
	bool useHDF5() const override {
		return useHDF5Flag;
	}

	/**
	 * To know if we should use the phase cut.
	 * \see IOptions.h
	 */
	bool usePhaseCut() const override {
		return usePhaseCutFlag;
	}

	/**
	 * Obtain the maximum value of impurities (He or Xe) to be used.
	 * \see IOptions.h
	 */
	int getMaxImpurity() const override {
		return maxImpurity;
	}

	/**
	 * Obtain the maximum value of deuterium to be used.
	 * \see IOptions.h
	 */
	int getMaxD() const override {
		return maxD;
	}

	/**
	 * Obtain the maximum value of tritium to be used.
	 * \see IOptions.h
	 */
	int getMaxT() const override {
		return maxT;
	}

	/**
	 * Obtain the maximum value of vacancies to be used.
	 * \see IOptions.h
	 */
	int getMaxV() const override {
		return maxV;
	}

	/**
	 * Obtain the maximum value of interstitials to be used.
	 * \see IOptions.h
	 */
	int getMaxI() const override {
		return maxI;
	}

	/**
	 * Obtain the number of grid points in the depth direction to be used.
	 * \see IOptions.h
	 */
	int getNX() const override {
		return nX;
	}

	/**
	 * Obtain the value of the step size in the depth direction to be used.
	 * \see IOptions.h
	 */
	double getXStepSize() const override {
		return xStepSize;
	}

	/**
	 * Obtain the number of grid points in the Y direction to be used.
	 * \see IOptions.h
	 */
	int getNY() const override {
		return nY;
	}

	/**
	 * Obtain the value of the step size in the Y direction to be used.
	 * \see IOptions.h
	 */
	double getYStepSize() const override {
		return yStepSize;
	}

	/**
	 * Obtain the number of grid points in the Z direction to be used.
	 * \see IOptions.h
	 */
	int getNZ() const override {
		return nZ;
	}

	/**
	 * Obtain the value of the step size in the Z direction to be used.
	 * \see IOptions.h
	 */
	double getZStepSize() const override {
		return zStepSize;
	}

	/**
	 * Obtain the boundary condition on a given side of the grid.
	 * \see IOptions.h
	 */
	int getLeftBoundary() const override {
		return leftBoundary;
	}
	int getRightBoundary() const override {
		return rightBoundary;
	}
	int getBottomBoundary() const override {
		return bottomBoundary;
	}
	int getTopBoundary() const override {
		return topBoundary;
	}
	int getFrontBoundary() const override {
		return frontBoundary;
	}
	int getBackBoundary() const override {
		return backBoundary;
	}

	/**
	 * Obtain the value of the depth above which the bursting is happening.
	 * \see IOptions.h
	 */
	double getBurstingDepth() const override {
		return burstingDepth;
	}

	/**
	 * Set the seed that should be used for initializing the random
	 * number generator.
	 * \see IOptions.h
	 */
	void setRNGSeed(unsigned int s) override {
		rngSeed = s;
		rngUseSeed = true;
	}

	/**
	 * Obtain the seed that should be used for initializing the random
	 * number generator.
	 * \see IOptions.h
	 */
	std::tuple<bool, unsigned int> getRNGSeed(void) const override {
		return std::make_tuple(rngUseSeed, rngSeed);
	}

	/**
	 * Determine if we should print the value used to seed the random
	 * number generator (regardless if it was given on the command line
	 * or generated dynamically).
	 * \see IOptions.h
	 */
	bool printRNGSeed(void) const override {
		return rngPrintSeed;
	}

	/**
	 * Obtain the minimum size for the re-solution.
	 * \see IOptions.h
	 */
	int getResoMinSize() const override {
		return resoMinSize;
	}

	/**
	 * Obtain the minimum size for the average radius computation.
	 * \see IOptions.h
	 */
	virtual util::Array<int, 4> getRadiusMinSizes() const override {
		return radiusMinSizes;
	}

	/**
	 * Obtain the value of the density of a bubble.
	 * \see IOptions.h
	 */
	double getDensity() const override {
		return density;
	}

	/**
	 * Obtain the value of the length of the flux pulse.
	 * \see IOptions.h
	 */
	virtual double getPulseTime() const override {
		return pulseTime;
	}

	/**
	 * Obtain the value of the proportion the flux pulse (on).
	 * \see IOptions.h
	 */
	virtual double getPulseProportion() const override {
		return pulseProportion;
	}

	/**
	 * Obtain the value of the lattice parameter.
	 * \see IOptions.h
	 */
	virtual double getLatticeParameter() const override {
		return latticeParameter;
	}

	/**
	 * Obtain the value of the impurity radius.
	 * \see IOptions.h
	 */
	virtual double getImpurityRadius() const override {
		return impurityRadius;
	}

	/**
	 * Obtain the value of the bias factor for interstitial.
	 * \see IOptions.h
	 */
	virtual double getBiasFactor() const override {
		return biasFactor;
	}

	/**
	 * Obtain the value of the factor between H and He radii.
	 * \see IOptions.h
	 */
	virtual double getHydrogenFactor() const override {
		return hydrogenFactor;
	}

	/**
	 * Obtain the value of the xenon diffusion coefficient.
	 * \see IOptions.h
	 */
	virtual double getXenonDiffusivity() const override {
		return xenonDiffusivity;
	}

	/**
	 * Obtain the value of the fission yield.
	 * \see IOptions.h
	 */
	virtual double getFissionYield() const override {
		return fissionYield;
	}

	/**
	 * Obtain the value of the migration threshold
	 * \see IOptions.h
	 */
	virtual double getMigrationThreshold() const override {
		return migrationThreshold;
	}

};
//end class Options

} /* namespace options */
} /* namespace xolotl */

#endif // OPTIONS_H
