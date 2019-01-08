#ifndef OPTIONS_H
#define OPTIONS_H

#include "IOptions.h"
#include <IOptionHandler.h>

namespace xolotlCore {

/**
 * Options realizes the IOptions interface.
 * All private members will be accessed through getters and setters.
 */
class Options: public IOptions {

protected:
	/**
	 * Map of options we support, keyed by option switch string.
	 */
	typedef std::map<std::string, IOptionHandler*> OptionsMap;
	OptionsMap optionsMap;

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
	 * The string of options that will be given to PETSc.
	 */
	std::string petscArgv;

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
	xolotlPerf::IHandlerRegistry::RegistryType perfRegistryType;

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
	void readParams(char* argv[]) override;

	/**
	 * Show our help message.
	 * \see IOptions.h
	 */
	void showHelp(std::ostream& os) const override;

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
	 * Set the name of the network file.
	 * \see IOptions.h
	 */
	void setNetworkFilename(const std::string& name) override {
		networkFilename = name;
	}

	/**
	 * Get the Argv for PETSc.
	 * \see IOptions.h
	 */
	std::string getPetscArgv() const override {
		return petscArgv;
	}

	/**
	 * Set the Argv for PETSc.
	 * \see IOptions.h
	 */
	void setPetscArgv(std::string argv) override {
		petscArgv = argv;
	}

	/**
	 * Should we use const temperature handlers?
	 * \see IOptions.h
	 */
	bool useConstTemperatureHandlers() const override {
		return constTempFlag;
	}

	/**
	 * Set the constTempFlag.
	 * \see IOptions.h
	 */
	void setConstTempFlag(bool flag) override {
		constTempFlag = flag;
	}

	/**
	 * Obtain the value of the constant temperature to be used.
	 * \see IOptions.h
	 */
	double getConstTemperature() const override {
		return constTemperature;
	}

	/**
	 * Set the constant temperature.
	 * \see IOptions.h
	 */
	void setConstTemperature(double temp) override {
		constTemperature = temp;
	}

	/**
	 * Should we use temperature profile handlers?
	 * \see IOptions.h
	 */
	bool useTemperatureProfileHandlers() const override {
		return tempProfileFlag;
	}

	/**
	 * Set the tempProfileFlag.
	 * \see IOptions.h
	 */
	void setTempProfileFlag(bool flag) override {
		tempProfileFlag = flag;
	}

	/**
	 * Obtain the name of the file containing the temperature profile data.
	 * \see IOptions.h
	 */
	std::string getTempProfileFilename() const override {
		return tempProfileFilename;
	}

	/**
	 * Set the name of the profile file to use.
	 * \see IOptions.h
	 */
	void setTempProfileFilename(const std::string& name) override {
		tempProfileFilename = name;
	}

	/**
	 * Should we use heat equation handlers?
	 * \see IOptions.h
	 */
	bool useHeatEquationHandlers() const override {
		return heatFlag;
	}

	/**
	 * Set the heatFlag.
	 * \see IOptions.h
	 */
	void setHeatFlag(bool flag) override {
		heatFlag = flag;
	}

	/**
	 * Obtain the value of the temperature to be used in the bulk.
	 * \see IOptions.h
	 */
	double getBulkTemperature() const override {
		return bulkTemperature;
	}

	/**
	 * Set the bulk temperature.
	 * \see IOptions.h
	 */
	void setBulkTemperature(double temp) override {
		bulkTemperature = temp;
	}

	/**
	 * Should we use the flux option?
	 * \see IOptions.h
	 */
	bool useFluxAmplitude() const override {
		return fluxFlag;
	}
	;

	/**
	 * Set the fluxFlag.
	 * \see IOptions.h
	 */
	void setFluxFlag(bool flag) override {
		fluxFlag = flag;
	}

	/**
	 * Obtain the value of the flux intensity to be used.
	 * \see IOptions.h
	 */
	double getFluxAmplitude() const override {
		return fluxAmplitude;
	}

	/**
	 * Set the value for the flux intensity to use.
	 * \see IOptions.h
	 */
	void setFluxAmplitude(double flux) override {
		fluxAmplitude = flux;
	}

	/**
	 * Should we use a time profile for the flux?
	 * \see IOptions.h
	 */
	bool useFluxTimeProfile() const override {
		return fluxProfileFlag;
	}

	/**
	 * Set the fluxProfileFlag.
	 * \see IOptions.h
	 */
	void setFluxProfileFlag(bool flag) override {
		fluxProfileFlag = flag;
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
	 * Set the name of the time profile file to use.
	 * \see IOptions.h
	 */
	void setFluxProfileName(const std::string& name) override {
		fluxProfileFilename = name;
	}

	/**
	 * Which type of performance handlers should we use?
	 * \see IOptions.h
	 */
	xolotlPerf::IHandlerRegistry::RegistryType getPerfHandlerType(void) const
			override {
		return perfRegistryType;
	}

	/**
	 * Set the type of performance handlers to use.
	 * \see IOptions.h
	 */
	void setPerfHandlerType(xolotlPerf::IHandlerRegistry::RegistryType rtype)
			override {
		perfRegistryType = rtype;
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
	 * Set the vizStandardHandlersFlag.
	 * \see IOptions.h
	 */
	void setVizStandardHandlers(bool flag) override {
		vizStandardHandlersFlag = flag;
	}

	/**
	 * Obtain the name of the material to be used for simulation.
	 * \see IOptions.h
	 */
	std::string getMaterial() const override {
		return materialName;
	}

	/**
	 * Set the name of the material to be used for the simulation.
	 * \see IOptions.h
	 */
	void setMaterial(const std::string& material) override {
		materialName = material;
	}

	/**
	 * Obtain the value of the concentration for the vacancies.
	 * \see IOptions.h
	 */
	double getInitialVConcentration() const override {
		return initialVConcentration;
	}

	/**
	 * Set the value of the concentration for the vacancies.
	 * \see IOptions.h
	 */
	void setInitialVConcentration(double conc) override {
		initialVConcentration = conc;
	}

	/**
	 * Obtain the number of dimensions for the simulation.
	 * \see IOptions.h
	 */
	int getDimensionNumber() const override {
		return dimensionNumber;
	}

	/**
	 * Set the number of dimensions for the simulation.
	 * \see IOptions.h
	 */
	void setDimensionNumber(int number) override {
		dimensionNumber = number;
	}

	/**
	 * Obtain the value of the void portion for the simulation.
	 * \see IOptions.h
	 */
	double getVoidPortion() const override {
		return voidPortion;
	}

	/**
	 * Set the value of the void portion for the surface to grow.
	 * \see IOptions.h
	 */
	void setVoidPortion(double portion) override {
		voidPortion = portion;
	}

	/**
	 * Should we use a regular grid on the x direction?
	 * \see IOptions.h
	 */
	bool useRegularXGrid() const override {
		return useRegularGridFlag;
	}

	/**
	 * Set the useRegularGridFlag.
	 * \see IOptions.h
	 */
	void setRegularXGrid(bool flag) override {
		useRegularGridFlag = flag;
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
	 * Set the physical process map.
	 *
	 * @param map The map
	 */
	void setProcesses(std::map<std::string, bool> map) override {
		processMap = map;
	}

	/**
	 * Obtain the string listing the wanted GB.
	 * \see IOptions.h
	 */
	std::string getGbString() const override {
		return gbList;
	}

	/**
	 * Set the string listing the wanted GB.
	 * \see IOptions.h
	 */
	void setGbString(const std::string& gbString) override {
		gbList = gbString;
	}

	/**
	 * Obtain the minimum size for the grouping.
	 * \see IOptions.h
	 */
	int getGroupingMin() const override {
		return groupingMin;
	}

	/**
	 * Set the minimum size for the grouping.
	 * \see IOptions.h
	 */
	void setGroupingMin(int size) override {
		groupingMin = size;
	}

	/**
	 * Obtain the first width for the grouping.
	 * \see IOptions.h
	 */
	int getGroupingWidthA() const override {
		return groupingWidthA;
	}

	/**
	 * Set the first width for the grouping.
	 * \see IOptions.h
	 */
	void setGroupingWidthA(int width) override {
		groupingWidthA = width;
	}

	/**
	 * Obtain the second width for the grouping.
	 * \see IOptions.h
	 */
	int getGroupingWidthB() const override {
		return groupingWidthB;
	}

	/**
	 * Set the second width for the grouping.
	 * \see IOptions.h
	 */
	void setGroupingWidthB(int width) override {
		groupingWidthB = width;
	}

	/**
	 * Obtain the value of the intensity of the sputtering yield to be used.
	 * \see IOptions.h
	 */
	double getSputteringYield() const override {
		return sputteringYield;
	}

	/**
	 * Set the value for the sputtering yield to use.
	 * \see IOptions.h
	 */
	void setSputteringYield(double yield) override {
		sputteringYield = yield;
	}

	/**
	 * To know if we should use the HDF5 file.
	 * \see IOptions.h
	 */
	bool useHDF5() const override {
		return useHDF5Flag;
	}

	/**
	 * Set the useHDF5Flag.
	 * \see IOptions.h
	 */
	void setHDF5Flag(bool flag) override {
		useHDF5Flag = flag;
	}

	/**
	 * To know if we should use the phase cut.
	 * \see IOptions.h
	 */
	bool usePhaseCut() const override {
		return usePhaseCutFlag;
	}

	/**
	 * Set the usePhaseCutFlag.
	 * \see IOptions.h
	 */
	void setPhaseCutFlag(bool flag) override {
		usePhaseCutFlag = flag;
	}

	/**
	 * Obtain the maximum value of impurities (He or Xe) to be used.
	 * \see IOptions.h
	 */
	int getMaxImpurity() const override {
		return maxImpurity;
	}

	/**
	 * Set the maximum value of impurities to use.
	 * \see IOptions.h
	 */
	void setMaxImpurity(int max) override {
		maxImpurity = max;
	}

	/**
	 * Obtain the maximum value of deuterium to be used.
	 * \see IOptions.h
	 */
	int getMaxD() const override {
		return maxD;
	}

	/**
	 * Set the maximum value of deuterium to use.
	 * \see IOptions.h
	 */
	void setMaxD(int max) override {
		maxD = max;
	}

	/**
	 * Obtain the maximum value of tritium to be used.
	 * \see IOptions.h
	 */
	int getMaxT() const override {
		return maxT;
	}

	/**
	 * Set the maximum value of tritium to use.
	 * \see IOptions.h
	 */
	void setMaxT(int max) override {
		maxT = max;
	}

	/**
	 * Obtain the maximum value of vacancies to be used.
	 * \see IOptions.h
	 */
	int getMaxV() const override {
		return maxV;
	}

	/**
	 * Set the maximum value of vacancies to use.
	 * \see IOptions.h
	 */
	void setMaxV(int max) override {
		maxV = max;
	}

	/**
	 * Obtain the maximum value of interstitials to be used.
	 * \see IOptions.h
	 */
	int getMaxI() const override {
		return maxI;
	}

	/**
	 * Set the maximum value of interstitials to use.
	 * \see IOptions.h
	 */
	void setMaxI(int max) override {
		maxI = max;
	}

	/**
	 * Obtain the number of grid points in the depth direction to be used.
	 * \see IOptions.h
	 */
	int getNX() const override {
		return nX;
	}

	/**
	 * Set the number of grid points in the depth direction to use.
	 * \see IOptions.h
	 */
	void setNX(int n) override {
		nX = n;
	}

	/**
	 * Obtain the value of the step size in the depth direction to be used.
	 * \see IOptions.h
	 */
	double getXStepSize() const override {
		return xStepSize;
	}

	/**
	 * Set the value for the step size in the depth direction to use.
	 * \see IOptions.h
	 */
	void setXStepSize(double stepSize) override {
		xStepSize = stepSize;
	}

	/**
	 * Obtain the number of grid points in the Y direction to be used.
	 * \see IOptions.h
	 */
	int getNY() const override {
		return nY;
	}

	/**
	 * Set the number of grid points in the Y direction to use.
	 * \see IOptions.h
	 */
	void setNY(int n) override {
		nY = n;
	}

	/**
	 * Obtain the value of the step size in the Y direction to be used.
	 * \see IOptions.h
	 */
	double getYStepSize() const override {
		return yStepSize;
	}

	/**
	 * Set the value for the step size in the Y direction to use.
	 * \see IOptions.h
	 */
	void setYStepSize(double stepSize) override {
		yStepSize = stepSize;
	}

	/**
	 * Obtain the number of grid points in the Z direction to be used.
	 * \see IOptions.h
	 */
	int getNZ() const override {
		return nZ;
	}

	/**
	 * Set the number of grid points in the Z direction to use.
	 * \see IOptions.h
	 */
	void setNZ(int n) override {
		nZ = n;
	}

	/**
	 * Obtain the value of the step size in the Z direction to be used.
	 * \see IOptions.h
	 */
	double getZStepSize() const override {
		return zStepSize;
	}

	/**
	 * Set the value for the step size in the Z direction to use.
	 * \see IOptions.h
	 */
	void setZStepSize(double stepSize) override {
		zStepSize = stepSize;
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
	 * Set the boundary condition on a given side of the grid.
	 * \see IOptions.h
	 */
	void setLeftBoundary(int n) override {
		leftBoundary = n;
	}
	void setRightBoundary(int n) override {
		rightBoundary = n;
	}
	void setBottomBoundary(int n) override {
		bottomBoundary = n;
	}
	void setTopBoundary(int n) override {
		topBoundary = n;
	}
	void setFrontBoundary(int n) override {
		frontBoundary = n;
	}
	void setBackBoundary(int n) override {
		backBoundary = n;
	}

	/**
	 * Obtain the value of the depth above which the bursting is happening.
	 * \see IOptions.h
	 */
	double getBurstingDepth() const override {
		return burstingDepth;
	}

	/**
	 * Set the value for the depth above which the bursting is happening.
	 * \see IOptions.h
	 */
	void setBurstingDepth(double depth) override {
		burstingDepth = depth;
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
	 * Specify whether each process should print the value it uses
	 * to seed the random number generator.
	 *
	 * @param b A bool indicating whether to print the RNG seed value.
	 */
	void setPrintRNGSeed(bool b) override {
		rngPrintSeed = b;
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

};
//end class Options

} /* namespace xolotlCore */

#endif // OPTIONS_H
