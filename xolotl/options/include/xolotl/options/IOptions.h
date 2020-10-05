#ifndef IOPTIONS_H
#define IOPTIONS_H

// Includes
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <xolotl/perf/xolotlPerf.h>
#include <xolotl/util/Array.h>
#include <xolotl/util/Filesystem.h>

namespace xolotl
{
namespace options
{
/**
 * IOptions describes the structure needed for the options in Xolotl.
 * All private members will be accessed through getters and setters.
 */
class IOptions
{
public:
	/**
	 * The destructor
	 */
	virtual ~IOptions()
	{
	}

	/**
	 * Read the parameters from the given file to set the different
	 * Xolotl options.
	 *
	 * @param argv Vector of argument strings
	 */
	virtual void
	readParams(int argc, char* argv[]) = 0;

	/**
	 * Should the program run after parsing the parameter file?
	 *
	 * @return true is the program should run
	 */
	virtual bool
	shouldRun() const = 0;

	/**
	 * Set the shouldRunFlag.
	 *
	 * @param flag The value for the shouldRunFlag
	 */
	virtual void
	setShouldRunFlag(bool flag) = 0;

	/**
	 * If program shouldn't run, what should its exit code be?
	 *
	 * @return the value of the exit code
	 */
	virtual int
	getExitCode() const = 0;

	/**
	 * Set the value for the exit code.
	 *
	 * @param code The value for exit code
	 */
	virtual void
	setExitCode(int code) = 0;

	/**
	 * Get the name of the network file.
	 *
	 * @return the name of the network file
	 */
	virtual std::string
	getNetworkFilename() const = 0;

	/**
	 * Get the name of the solver to use
	 *
	 * @return the name of the solver
	 */
	virtual std::string
	getSolverName() const = 0;

	/**
	 * Get the Arguments for PETSc.
	 *
	 * @return arg
	 */
	virtual std::string
	getPetscArg() const = 0;

	/**
	 * Obtain the name of the temperature handler to be used
	 *
	 * @return The name of the temperature handler
	 */
	virtual std::string
	getTempHandlerName() const = 0;

	/**
	 * Obtain the temperature parameters for the temperature handler
	 */
	virtual double
	getTempParam(std::size_t i = 0) const = 0;

	/**
	 * Obtain the name of the file containing the temperature profile data.
	 *
	 * @return The name of the file
	 */
	virtual std::string
	getTempProfileFilename() const = 0;

	/**
	 * Should we use the flux amplitude option?
	 * If false, it will not be used.
	 *
	 * @return true if the flux amplitude option was present in
	 * the parameter file, false if it was not
	 */
	virtual bool
	useFluxAmplitude() const = 0;

	/**
	 * Obtain the value of the intensity of the flux to be used.
	 *
	 * @return The value of the flux amplitude
	 */
	virtual double
	getFluxAmplitude() const = 0;

	/**
	 * Should we use a time profile for the helium flux?
	 *
	 * @return True is a time profile file is given for the helium flux
	 */
	virtual bool
	useFluxTimeProfile() const = 0;

	/**
	 * Obtain the name of the file containing the time profile data for the
	 * helium flux.
	 *
	 * @return The name of the file
	 */
	virtual std::string
	getFluxTimeProfileFilePath() const = 0;

	/**
	 * Which type of performance handlers should we use?
	 *
	 * @return The type of performance handler registry to use
	 */
	virtual perf::IHandlerRegistry::RegistryType
	getPerfHandlerType(void) const = 0;

	/**
	 * Obtain the name of the visualization handler to be used
	 *
	 * @return The name of the viz handler
	 */
	virtual std::string
	getVizHandlerName() const = 0;

	/**
	 * Obtain the name of the material to be used for the simulation.
	 *
	 * @return The name of the material
	 */
	virtual std::string
	getMaterial() const = 0;

	/**
	 * Obtain the value of the void portion for the simulation.
	 * @return The portion.
	 */
	virtual double
	getVoidPortion() const = 0;

	/**
	 * Obtain the value of the concentration for the vacancies.
	 *
	 * @return The concentration value
	 */
	virtual double
	getInitialVConcentration() const = 0;

	/**
	 * Obtain the value of the electronic stopping power.
	 *
	 * @return Zeta
	 */
	virtual double
	getZeta() const = 0;

	/**
	 * Obtain the number of dimensions for the simulation.
	 *
	 * @return The number of dimensions
	 */
	virtual int
	getDimensionNumber() const = 0;

	/**
	 * Should we use a regular grid on the x direction?
	 * @return true if program should use a regular grid,
	 * false if not
	 */
	virtual bool
	useRegularXGrid() const = 0;

	/**
	 * Should we use a Chebyshev grid on the x direction?
	 * @return true if program should use a Chebyshev grid,
	 * false if not
	 */
	virtual bool
	useChebyshevGrid() const = 0;

	/**
	 * Should we read in the grid on the x direction?
	 * @return true if program should read in the grid,
	 * false if not
	 */
	virtual bool
	useReadInGrid() const = 0;

	/**
	 * Get the name of the grid file.
	 *
	 * @return The name of the grid file
	 */
	virtual std::string
	getGridFilename() const = 0;

	/**
	 * Obtain the physical process map.
	 *
	 * @return The map
	 */
	virtual const std::map<std::string, bool>&
	getProcesses() const = 0;

	/**
	 * Obtain the string listing the wanted GB.
	 *
	 * @return The string of GB
	 */
	virtual std::string
	getGbString() const = 0;

	/**
	 * Obtain the minimum size for the grouping.
	 *
	 * @return The size
	 */
	virtual int
	getGroupingMin() const = 0;

	/**
	 * Obtain the first width for the grouping.
	 *
	 * @return The width
	 */
	virtual int
	getGroupingWidthA() const = 0;

	/**
	 * Obtain the second width for the grouping.
	 *
	 * @return The width
	 */
	virtual int
	getGroupingWidthB() const = 0;

	/**
	 * Obtain the value of the intensity of the sputtering yield to be used.
	 *
	 * @return The value of the sputtering yield
	 */
	virtual double
	getSputteringYield() const = 0;

	/**
	 * To know if we should use the HDF5 file.
	 *
	 * @return useHDF5Flag
	 */
	virtual bool
	useHDF5() const = 0;

	/**
	 * To know if we should use the phase cut.
	 *
	 * @return usePhaseCutFlag
	 */
	virtual bool
	usePhaseCut() const = 0;

	/**
	 * Obtain the maximum value of impurities (He or Xe) to be used.
	 *
	 * @return The maximum value
	 */
	virtual int
	getMaxImpurity() const = 0;

	/**
	 * Obtain the maximum value of deuterium to be used.
	 *
	 * @return The maximum value
	 */
	virtual int
	getMaxD() const = 0;

	/**
	 * Obtain the maximum value of tritium to be used.
	 *
	 * @return The maximum value
	 */
	virtual int
	getMaxT() const = 0;

	/**
	 * Obtain the maximum value of vacancies to be used.
	 *
	 * @return The maximum value
	 */
	virtual int
	getMaxV() const = 0;

	/**
	 * Obtain the maximum value of interstitials to be used.
	 *
	 * @return The maximum value
	 */
	virtual int
	getMaxI() const = 0;

	/**
	 * Obtain the number of grid points in the depth direction to be used.
	 *
	 * @return The number of grid points
	 */
	virtual int
	getNX() const = 0;

	/**
	 * Obtain the value of the step size in the depth direction to be used.
	 *
	 * @return The value of the step size
	 */
	virtual double
	getXStepSize() const = 0;

	/**
	 * Obtain the number of grid points in the Y direction to be used.
	 *
	 * @return The number of grid points
	 */
	virtual int
	getNY() const = 0;

	/**
	 * Obtain the value of the step size in the Y direction to be used.
	 *
	 * @return The value of the step size
	 */
	virtual double
	getYStepSize() const = 0;

	/**
	 * Obtain the number of grid points in the Z direction to be used.
	 *
	 * @return The number of grid points
	 */
	virtual int
	getNZ() const = 0;

	/**
	 * Obtain the value of the step size in the Z direction to be used.
	 *
	 * @return The value of the step size
	 */
	virtual double
	getZStepSize() const = 0;

	/**
	 * Obtain the boundary condition on a given side of the grid.
	 *
	 * @return The boundary condition
	 */
	virtual int
	getLeftBoundary() const = 0;
	virtual int
	getRightBoundary() const = 0;
	virtual int
	getBottomBoundary() const = 0;
	virtual int
	getTopBoundary() const = 0;
	virtual int
	getFrontBoundary() const = 0;
	virtual int
	getBackBoundary() const = 0;

	/**
	 * Obtain the value of the depth above which the bursting is happening.
	 *
	 * @return The depth
	 */
	virtual double
	getBurstingDepth() const = 0;

	/**
	 * Obtain the value of the factor in the bursting probability.
	 *
	 * @return The factor
	 */
	virtual double
	getBurstingFactor() const = 0;

	/**
	 * Set the seed that should be used for initializing the random
	 * number generator.
	 *
	 * @param s The value to use to seed the RNG.
	 */
	virtual void
	setRNGSeed(unsigned int s) = 0;

	/**
	 * Obtain the seed that should be used for initializing the random
	 * number generator.
	 *
	 * @return A (bool, uint) pair.  The bool tells whether to use the int
	 * to seed the random number generator.
	 */
	virtual std::tuple<bool, unsigned int>
	getRNGSeed(void) const = 0;

	/**
	 * Determine if we should print the value used to seed the random
	 * number generator (regardless if it was given on the command line
	 * or generated dynamically).
	 *
	 * @return True iff we should print the RNG seed value from each process.
	 */
	virtual bool
	printRNGSeed(void) const = 0;

	/**
	 * Obtain the minimum size for the average radius computation.
	 *
	 * @return The size
	 */
	virtual std::vector<size_t>
	getRadiusMinSizes() const = 0;

	/**
	 * Obtain the value of the density of a bubble.
	 *
	 * @return The density
	 */
	virtual double
	getDensity() const = 0;

	/**
	 * Obtain the value of the length of the flux pulse.
	 *
	 * @return The value for the pulse
	 */
	virtual double
	getPulseTime() const = 0;

	/**
	 * Obtain the value of the proportion the flux pulse (on).
	 *
	 * @return The value for the proportion
	 */
	virtual double
	getPulseProportion() const = 0;

	/**
	 * Obtain the value of the lattice parameter.
	 *
	 * @return The lattice parameter in nm
	 */
	virtual double
	getLatticeParameter() const = 0;

	/**
	 * Obtain the value of the impurity radius.
	 *
	 * @return The impurity radius in nm
	 */
	virtual double
	getImpurityRadius() const = 0;

	/**
	 * Obtain the value of the bias factor for interstitial.
	 *
	 * @return The bias
	 */
	virtual double
	getBiasFactor() const = 0;

	/**
	 * Obtain the value of the factor between H and He radii.
	 *
	 * @return The hydrogen factor
	 */
	virtual double
	getHydrogenFactor() const = 0;

	/**
	 * Obtain the value of the xenon diffusion coefficient.
	 *
	 * @return The diffusivity in nm2 s-1
	 */
	virtual double
	getXenonDiffusivity() const = 0;

	/**
	 * Obtain the value of the fission yield.
	 *
	 * @return The number of xenon per fission
	 */
	virtual double
	getFissionYield() const = 0;

	/**
	 * Obtain the value of the HeV ratio.
	 *
	 * @return The ratio
	 */
	virtual double
	getHeVRatio() const = 0;

	/**
	 * Obtain the value of the migration energy threshold for effective
	 * diffusivity.
	 *
	 * @return The value of the threshold
	 */
	virtual double
	getMigrationThreshold() const = 0;

	/**
	 * Get the path to the custom flux profile.
	 *
	 * @return The path to the file
	 */
	virtual std::string
	getFluxDepthProfileFilePath() const = 0;
};
// end class IOptions

} /* namespace options */
} /* namespace xolotl */

#endif
