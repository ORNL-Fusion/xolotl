#ifndef IOPTIONS_H
#define IOPTIONS_H

// Includes
#include <iostream>
#include <string>
#include <map>
#include <xolotlPerf.h>
#include <NDArray.h>

namespace xolotlCore {

/**
 * IOptions describes the structure needed for the options in Xolotl.
 * All private members will be accessed through getters and setters.
 */
class IOptions {

public:

	/**
	 * The destructor
	 */
	virtual ~IOptions() {
	}

	/**
	 * Read the parameters from the given file to set the different
	 * Xolotl options.
	 *
	 * @param argv Vector of argument strings
	 */
	virtual void readParams(int argc, char* argv[]) = 0;

	/**
	 * Should the program run after parsing the parameter file?
	 *
	 * @return true is the program should run
	 */
	virtual bool shouldRun() const = 0;

	/**
	 * Set the shouldRunFlag.
	 *
	 * @param flag The value for the shouldRunFlag
	 */
	virtual void setShouldRunFlag(bool flag) = 0;

	/**
	 * If program shouldn't run, what should its exit code be?
	 *
	 * @return the value of the exit code
	 */
	virtual int getExitCode() const = 0;

	/**
	 * Set the value for the exit code.
	 *
	 * @param code The value for exit code
	 */
	virtual void setExitCode(int code) = 0;

	/**
	 * Get the name of the network file.
	 *
	 * @return the name of the network file
	 */
	virtual std::string getNetworkFilename() const = 0;

	/**
	 * Set the name of the network file.
	 *
	 * @param name Name for the network file
	 */
	virtual void setNetworkFilename(const std::string& name) = 0;

	/**
	 * Get the Argc for PETSc.
	 *
	 * @return argc
	 */
	virtual int getPetscArgc() const = 0;

	/**
	 * Set the Argc for PETSc.
	 *
	 * @param argc The number of options for PETSc
	 */
	virtual void setPetscArgc(int argc) = 0;

	/**
	 * Get the Argv for PETSc.
	 *
	 * @return argv
	 */
	virtual char** getPetscArgv() const = 0;

	/**
	 * Set the Argv for PETSc.
	 *
	 * @param argv The pointer to the options for PETSc
	 */
	virtual void setPetscArgv(char** argv) = 0;

	/**
	 * Should we use const temperature handlers?
	 *
	 * @return true if Xolotl must use a constant temperature
	 */
	virtual bool useConstTemperatureHandlers() const = 0;

	/**
	 * Set the constTempFlag.
	 *
	 * @param flag The value for the constTempFlag
	 */
	virtual void setConstTempFlag(bool flag) = 0;

	/**
	 * Obtain the value of the constant temperature to be used.
	 *
	 * @return The value for the temperature
	 */
	virtual double getConstTemperature() const = 0;

	/**
	 * Set the constant temperature.
	 *
	 * @param temp The value for the constant temperature
	 */
	virtual void setConstTemperature(double temp) = 0;

	/**
	 * Should we use temperature profile handlers?
	 *
	 * @return true if Xolotl must use a temperature profile
	 */
	virtual bool useTemperatureProfileHandlers() const = 0;

	/**
	 * Set the tempProfileFlag.
	 *
	 * @param flag The value for the tempProfileFlag
	 */
	virtual void setTempProfileFlag(bool flag) = 0;

	/**
	 * Obtain the name of the file containing the temperature profile data.
	 *
	 * @return The name of the file
	 */
	virtual std::string getTempProfileFilename() const = 0;

	/**
	 * Set the name of the profile file to use.
	 *
	 * @param name The name of the file
	 */
	virtual void setTempProfileFilename(const std::string& name) = 0;

	/**
	 * Should we use heat equation handlers?
	 *
	 * @return true if Xolotl must use the heat equation
	 */
	virtual bool useHeatEquationHandlers() const = 0;

	/**
	 * Set the heatFlag.
	 *
	 * @param flag The value for the heatFlag
	 */
	virtual void setHeatFlag(bool flag) = 0;

	/**
	 * Obtain the value of the temperature to be used in the bulk.
	 *
	 * @return The value for the temperature
	 */
	virtual double getBulkTemperature() const = 0;

	/**
	 * Set the bulk temperature.
	 *
	 * @param temp The value for the bulk temperature
	 */
	virtual void setBulkTemperature(double temp) = 0;

	/**
	 * Should we use the flux amplitude option?
	 * If false, it will not be used.
	 *
	 * @return true if the flux amplitude option was present in
	 * the parameter file, false if it was not
	 */
	virtual bool useFluxAmplitude() const = 0;

	/**
	 * Set the fluxFlag.
	 *
	 * @param flag The value for the fluxFlag
	 */
	virtual void setFluxFlag(bool flag) = 0;

	/**
	 * Obtain the value of the intensity of the flux to be used.
	 *
	 * @return The value of the flux amplitude
	 */
	virtual double getFluxAmplitude() const = 0;

	/**
	 * Set the value for the flux intensity to use.
	 *
	 * @param flux The value for the flux amplitude
	 */
	virtual void setFluxAmplitude(double flux) = 0;

	/**
	 * Should we use a time profile for the helium flux?
	 *
	 * @return True is a time profile file is given for the helium flux
	 */
	virtual bool useFluxTimeProfile() const = 0;

	/**
	 * Set the fluxProfileFlag.
	 *
	 * @param flag The value for the flag
	 */
	virtual void setFluxProfileFlag(bool flag) = 0;

	/**
	 * Obtain the name of the file containing the time profile data for the
	 * helium flux.
	 *
	 * @return The name of the file
	 */
	virtual std::string getFluxProfileName() const = 0;

	/**
	 * Set the name of the time profile file to use.
	 *
	 * @param name The name of the file
	 */
	virtual void setFluxProfileName(const std::string& name) = 0;

	/**
	 * Which type of performance handlers should we use?
	 *
	 * @return The type of performance handler registry to use
	 */
	virtual xolotlPerf::IHandlerRegistry::RegistryType getPerfHandlerType(
			void) const = 0;

	/**
	 * Set the type of performance handlers to use.
	 *
	 * @param rtype The type of performance handler registry to use
	 */
	virtual void setPerfHandlerType(
			xolotlPerf::IHandlerRegistry::RegistryType rtype) = 0;

	/**
	 * Should we use the "standard" set of handlers for the visualization?
	 * If false, use dummy (stub) handlers.
	 *
	 * @return true if program should use standard handlers, false if
	 * should use dummy handlers
	 */
	virtual bool useVizStandardHandlers() const = 0;

	/**
	 * Set the vizStandardHandlersFlag.
	 *
	 * @param flag The value for the vizStandardHandlersFlag
	 */
	virtual void setVizStandardHandlers(bool flag) = 0;

	/**
	 * Obtain the name of the material to be used for the simulation.
	 *
	 * @return The name of the material
	 */
	virtual std::string getMaterial() const = 0;

	/**
	 * Set the name of the material to be used for the simulation.
	 *
	 * @param material The name of the material
	 */
	virtual void setMaterial(const std::string& material) = 0;

	/**
	 * Obtain the value of the void portion for the simulation.
	 * @return The portion.
	 */
	virtual double getVoidPortion() const = 0;

	/**
	 * Set the value of the void portion for the surface to grow.
	 * @param portion The value for the portion.
	 */
	virtual void setVoidPortion(double portion) = 0;

	/**
	 * Obtain the value of the concentration for the vacancies.
	 *
	 * @return The concentration value
	 */
	virtual double getInitialVConcentration() const = 0;

	/**
	 * Set the value of the concentration for the vacancies.
	 *
	 * @param conc The value for the concentration
	 */
	virtual void setInitialVConcentration(double conc) = 0;

	/**
	 * Obtain the value of the electronic stopping power.
	 *
	 * @return Zeta
	 */
	virtual double getZeta() const = 0;

	/**
	 * Set the value of the electronic stopping power.
	 *
	 * @param zeta The value for the electronic stopping power
	 */
	virtual void setZeta(double zeta) = 0;

	/**
	 * Obtain the number of dimensions for the simulation.
	 *
	 * @return The number of dimensions
	 */
	virtual int getDimensionNumber() const = 0;

	/**
	 * Set the number of dimensions for the simulation.
	 *
	 * @param number The number of dimensions
	 */
	virtual void setDimensionNumber(int number) = 0;

	/**
	 * Should we use a regular grid on the x direction?
	 * @return true if program should use a regular grid,
	 * false if not
	 */
	virtual bool useRegularXGrid() const = 0;

	/**
	 * Set the useRegularGridFlag.
	 * @param flag The value for the useRegularGridFlag.
	 */
	virtual void setRegularXGrid(bool flag) = 0;

	/**
	 * Should we use a Chebyshev grid on the x direction?
	 * @return true if program should use a Chebyshev grid,
	 * false if not
	 */
	virtual bool useChebyshevGrid() const = 0;

	/**
	 * Set the useChebyshevGridFlag.
	 * @param flag The value for the useChebyshevGridFlag.
	 */
	virtual void setChebyshevGrid(bool flag) = 0;

	/**
	 * Should we read in the grid on the x direction?
	 * @return true if program should read in the grid,
	 * false if not
	 */
	virtual bool useReadInGrid() const = 0;

	/**
	 * Set the readInGridFlag.
	 * @param flag The value for the readInGridFlag.
	 */
	virtual void setReadInGrid(bool flag) = 0;

	/**
	 * Get the name of the grid file.
	 *
	 * @return The name of the grid file
	 */
	virtual std::string getGridFilename() const = 0;

	/**
	 * Set the name of the grid file.
	 *
	 * @param name Name for the grid file
	 */
	virtual void setGridFilename(const std::string& name) = 0;

	/**
	 * Obtain the physical process map.
	 *
	 * @return The map
	 */
	virtual std::map<std::string, bool> getProcesses() const = 0;

	/**
	 * Set the physical process map.
	 *
	 * @param map The map
	 */
	virtual void setProcesses(std::map<std::string, bool> map) = 0;

	/**
	 * Obtain the string listing the wanted GB.
	 *
	 * @return The string of GB
	 */
	virtual std::string getGbString() const = 0;

	/**
	 * Set the string listing the wanted GB.
	 *
	 * @param gbString The string of GB
	 */
	virtual void setGbString(const std::string& gbString) = 0;

	/**
	 * Obtain the minimum size for the grouping.
	 *
	 * @return The size
	 */
	virtual int getGroupingMin() const = 0;

	/**
	 * Set the minimum size for the grouping.
	 *
	 * @param size The size
	 */
	virtual void setGroupingMin(int size) = 0;

	/**
	 * Obtain the first width for the grouping.
	 *
	 * @return The width
	 */
	virtual int getGroupingWidthA() const = 0;

	/**
	 * Set the first width for the grouping.
	 *
	 * @param width The width
	 */
	virtual void setGroupingWidthA(int width) = 0;

	/**
	 * Obtain the second width for the grouping.
	 *
	 * @return The width
	 */
	virtual int getGroupingWidthB() const = 0;

	/**
	 * Set the second width for the grouping.
	 *
	 * @param width The width
	 */
	virtual void setGroupingWidthB(int width) = 0;

	/**
	 * Obtain the value of the intensity of the sputtering yield to be used.
	 *
	 * @return The value of the sputtering yield
	 */
	virtual double getSputteringYield() const = 0;

	/**
	 * Set the value for the sputtering yield to use.
	 *
	 * @param yield The value for the sputtering yield
	 */
	virtual void setSputteringYield(double yield) = 0;

	/**
	 * To know if we should use the HDF5 file.
	 *
	 * @return useHDF5Flag
	 */
	virtual bool useHDF5() const = 0;

	/**
	 * Set the useHDF5Flag.
	 *
	 * @param flag The value for the useHDF5Flag
	 */
	virtual void setHDF5Flag(bool flag) = 0;

	/**
	 * To know if we should use the phase cut.
	 *
	 * @return usePhaseCutFlag
	 */
	virtual bool usePhaseCut() const = 0;

	/**
	 * Set the usePhaseCutFlag.
	 *
	 * @param flag The value for the usePhaseCutFlag
	 */
	virtual void setPhaseCutFlag(bool flag) = 0;

	/**
	 * Obtain the maximum value of impurities (He or Xe) to be used.
	 *
	 * @return The maximum value
	 */
	virtual int getMaxImpurity() const = 0;

	/**
	 * Set the maximum value of impurities to use.
	 *
	 * @param max The maximum
	 */
	virtual void setMaxImpurity(int max) = 0;

	/**
	 * Obtain the maximum value of deuterium to be used.
	 *
	 * @return The maximum value
	 */
	virtual int getMaxD() const = 0;

	/**
	 * Set the maximum value of deuterium to use.
	 *
	 * @param max The maximum
	 */
	virtual void setMaxD(int max) = 0;

	/**
	 * Obtain the maximum value of tritium to be used.
	 *
	 * @return The maximum value
	 */
	virtual int getMaxT() const = 0;

	/**
	 * Set the maximum value of tritium to use.
	 *
	 * @param max The maximum
	 */
	virtual void setMaxT(int max) = 0;

	/**
	 * Obtain the maximum value of vacancies to be used.
	 *
	 * @return The maximum value
	 */
	virtual int getMaxV() const = 0;

	/**
	 * Set the maximum value of vacancies to use.
	 *
	 * @param max The maximum
	 */
	virtual void setMaxV(int max) = 0;

	/**
	 * Obtain the maximum value of interstitials to be used.
	 *
	 * @return The maximum value
	 */
	virtual int getMaxI() const = 0;

	/**
	 * Set the maximum value of interstitials to use.
	 *
	 * @param max The maximum
	 */
	virtual void setMaxI(int max) = 0;

	/**
	 * Obtain the number of grid points in the depth direction to be used.
	 *
	 * @return The number of grid points
	 */
	virtual int getNX() const = 0;

	/**
	 * Set the number of grid points in the depth direction to use.
	 *
	 * @param n The number
	 */
	virtual void setNX(int n) = 0;

	/**
	 * Obtain the value of the step size in the depth direction to be used.
	 *
	 * @return The value of the step size
	 */
	virtual double getXStepSize() const = 0;

	/**
	 * Set the value for the step size in the depth direction to use.
	 *
	 * @param stepSize The value for the step size
	 */
	virtual void setXStepSize(double stepSize) = 0;

	/**
	 * Obtain the number of grid points in the Y direction to be used.
	 *
	 * @return The number of grid points
	 */
	virtual int getNY() const = 0;

	/**
	 * Set the number of grid points in the Y direction to use.
	 *
	 * @param n The number
	 */
	virtual void setNY(int n) = 0;

	/**
	 * Obtain the value of the step size in the Y direction to be used.
	 *
	 * @return The value of the step size
	 */
	virtual double getYStepSize() const = 0;

	/**
	 * Set the value for the step size in the Y direction to use.
	 *
	 * @param stepSize The value for the step size
	 */
	virtual void setYStepSize(double stepSize) = 0;

	/**
	 * Obtain the number of grid points in the Z direction to be used.
	 *
	 * @return The number of grid points
	 */
	virtual int getNZ() const = 0;

	/**
	 * Set the number of grid points in the Z direction to use.
	 *
	 * @param n The number
	 */
	virtual void setNZ(int n) = 0;

	/**
	 * Obtain the value of the step size in the Z direction to be used.
	 *
	 * @return The value of the step size
	 */
	virtual double getZStepSize() const = 0;

	/**
	 * Set the value for the step size in the Z direction to use.
	 *
	 * @param stepSize The value for the step size
	 */
	virtual void setZStepSize(double stepSize) = 0;

	/**
	 * Obtain the boundary condition on a given side of the grid.
	 *
	 * @return The boundary condition
	 */
	virtual int getLeftBoundary() const = 0;
	virtual int getRightBoundary() const = 0;
	virtual int getBottomBoundary() const = 0;
	virtual int getTopBoundary() const = 0;
	virtual int getFrontBoundary() const = 0;
	virtual int getBackBoundary() const = 0;

	/**
	 * Set the boundary condition on a given side of the grid.
	 *
	 * @param n The condition
	 */
	virtual void setLeftBoundary(int n) = 0;
	virtual void setRightBoundary(int n) = 0;
	virtual void setBottomBoundary(int n) = 0;
	virtual void setTopBoundary(int n) = 0;
	virtual void setFrontBoundary(int n) = 0;
	virtual void setBackBoundary(int n) = 0;

	/**
	 * Obtain the value of the depth above which the bursting is happening.
	 *
	 * @return The depth
	 */
	virtual double getBurstingDepth() const = 0;

	/**
	 * Set the value for the depth above which the bursting is happening.
	 *
	 * @param depth The depth
	 */
	virtual void setBurstingDepth(double depth) = 0;

	/**
	 * Set the seed that should be used for initializing the random
	 * number generator.
	 *
	 * @param s The value to use to seed the RNG.
	 */
	virtual void setRNGSeed(unsigned int s) = 0;

	/**
	 * Obtain the seed that should be used for initializing the random
	 * number generator.
	 *
	 * @return A (bool, uint) pair.  The bool tells whether to use the int
	 * to seed the random number generator.
	 */
	virtual std::tuple<bool, unsigned int> getRNGSeed(void) const = 0;

	/**
	 * Specify whether each process should print the value it uses
	 * to seed the random number generator.
	 *
	 * @param b A bool indicating whether to print the RNG seed value.
	 */
	virtual void setPrintRNGSeed(bool b) = 0;

	/**
	 * Determine if we should print the value used to seed the random
	 * number generator (regardless if it was given on the command line
	 * or generated dynamically).
	 *
	 * @return True iff we should print the RNG seed value from each process.
	 */
	virtual bool printRNGSeed(void) const = 0;

	/**
	 * Obtain the minimum size for the re-solution.
	 *
	 * @return The size
	 */
	virtual int getResoMinSize() const = 0;

	/**
	 * Set the minimum size for the re-solution.
	 *
	 * @param size The minimum size
	 */
	virtual void setResoMinSize(int size) = 0;

	/**
	 * Obtain the minimum size for the average radius computation.
	 *
	 * @return The size
	 */
	virtual Array<int, 4> getRadiusMinSizes() const = 0;

	/**
	 * Set the minimum size for the average radius computation.
	 *
	 * @param sizes The minimum sizes
	 */
	virtual void setRadiusMinSizes(Array<int, 4> sizes) = 0;

	/**
	 * Obtain the value of the density of a bubble.
	 *
	 * @return The density
	 */
	virtual double getDensity() const = 0;

	/**
	 * Set the value for the density of a bubble.
	 *
	 * @param rho The density
	 */
	virtual void setDensity(double rho) = 0;

};
//end class IOptions

} /* namespace xolotlCore */

#endif
