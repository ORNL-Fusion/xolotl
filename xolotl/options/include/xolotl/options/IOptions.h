#pragma once

// Includes
#include <array>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include <xolotl/config.h>
#include <xolotl/util/Array.h>
#include <xolotl/util/Filesystem.h>

namespace xolotl
{
namespace interface
{
class MultiXolotl;
}
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

	virtual void
	printAll(std::ostream& os) const = 0;

	/**
	 * Create copy of all options
	 */
	virtual std::shared_ptr<IOptions>
	makeCopy() const = 0;

	/**
	 * Read the parameters from the given file to set the different
	 * Xolotl options.
	 *
	 * @param argc, argv Argument strings
	 */
	virtual void
	readParams(int argc, const char* argv[]) = 0;

	/**
	 * Get the name of the checkpoint file.
	 *
	 * @return the name of the checkpoint file
	 */
	virtual std::string
	getCheckpointFilePath() const = 0;

	/**
	 * Get the name of the restart file (could be different from checkpoint).
	 *
	 * @return the name of the restart file
	 */
	virtual std::string
	getRestartFilePath() const = 0;

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
	 * Obtain the temperature grid power.
	 */
	virtual double
	getTempGridPower() const = 0;

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
	 * Should we use a time profile for the flux?
	 *
	 * @return True is a time profile file is given for the flux
	 */
	virtual bool
	useFluxTimeProfile() const = 0;

	/**
	 * Obtain the name of the file containing the time profile data for the
	 * flux.
	 *
	 * @return The name of the file
	 */
	virtual std::string
	getFluxTimeProfileFilePath() const = 0;

	/**
	 * Obtain the name of the perfomance handler to be used
	 *
	 * @return The name of the perf handler
	 */
	virtual std::string
	getPerfHandlerName() const = 0;

	/**
	 * Should we write the performance report to a YAML file?
	 *
	 * @return true to enable YAML output
	 */
	virtual bool
	usePerfOutputYAML() const = 0;

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
	 * Obtain the interface location in nm.
	 * @return The location.
	 */
	virtual double
	getInterfaceLocation() const = 0;

	/**
	 * Obtain the value of the concentration for the vacancies.
	 *
	 * @return The concentration value
	 */
	virtual std::string
	getInitialConcentration() const = 0;

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
	 * Obtain the name of the grid type to be used
	 *
	 * @return The name of the grid type
	 */
	virtual std::string
	getGridTypeName() const = 0;

	/**
	 * Obtain the grid parameters
	 */
	virtual double
	getGridParam(std::size_t i = 0) const = 0;

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
	 * Add enabled entry to process map
	 */
	virtual void
	addProcess(const std::string& processKey) = 0;

	/**
	 * Determine if the network is being distributed to subnetworks
	 */
	virtual bool
	useSubnetworks() const = 0;

	/**
	 * Obtain the initial coupling time step
	 *
	 * @return The initial value
	 */
	virtual double
	getInitialTimeStep() const = 0;

	/**
	 * Obtain the maximum coupling time step
	 *
	 * @return The maximum value
	 */
	virtual double
	getMaxTimeStep() const = 0;

	/**
	 * Obtain the coupling time step growth factor
	 *
	 * @return The growth factor
	 */
	virtual double
	getTimeStepGrowthFactor() const = 0;

	/**
	 * Obtain the coupling start time
	 *
	 * @return The start time
	 */
	virtual double
	getStartTime() const = 0;

	/**
	 * Obtain the coupling end time
	 *
	 * @return The end time
	 */
	virtual double
	getEndTime() const = 0;

	/**
	 * Obtain the max number of coupling time steps
	 *
	 * @return The max number of steps
	 */
	virtual IdType
	getNumberOfTimeSteps() const = 0;

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
	 * Obtain the list of network parameters
	 *
	 * @return vector of IDs
	 */
	virtual const std::vector<IdType>&
	getNetworkParameters() const = 0;

	/**
	 * Replace all network parameters with given list
	 *
	 * @param params List of network parameters
	 */
	virtual void
	setNetworkParameters(const std::vector<IdType>& params) = 0;

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
	 * Obtain the maximum value of pure vacancies to be used.
	 *
	 * @return The maximum value
	 */
	virtual int
	getMaxPureV() const = 0;

	/**
	 * Obtain the maximum value of interstitials to be used.
	 *
	 * @return The maximum value
	 */
	virtual int
	getMaxI() const = 0;

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
	 * Obtain the string listing the wanted BC in the X direction.
	 *
	 * @return The BC type
	 */
	virtual std::string
	getBCString() const = 0;

	/**
	 * Obtain the value of the portion of heat lost to the bulk.
	 *
	 * @return The portion
	 */
	virtual double
	getHeatLossPortion() const = 0;

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

	/**
	 * Get the path to the reaction rate file.
	 *
	 * @return The path to the file
	 */
	virtual std::string
	getReactionFilePath() const = 0;

	/**
	 * Obtain the value of the basal portion.
	 *
	 * @return Qb
	 */
	virtual double
	getBasalPortion() const = 0;

	/**
	 * Obtain the transition size.
	 *
	 * @return The size
	 */
	virtual int
	getTransitionSize() const = 0;

	/**
	 * Obtain the value of dose at which the cascade overlap takes effect.
	 *
	 * @return The dose
	 */
	virtual double
	getCascadeDose() const = 0;

	/**
	 * Obtain the value of remaining cascade efficiency.
	 *
	 * @return The efficiency
	 */
	virtual double
	getCascadeEfficiency() const = 0;

protected:
	friend class ::xolotl::interface::MultiXolotl;

	/**
	 * Set checkpoint file (should only be used by MultiXolotl)
	 */
	virtual void
	setCheckpointFilePath(const std::string& path) = 0;

	/**
	 * Set restart file
	 *
	 * This should only be used by MultiXolotl. For a single instance, the file
	 * name is read from the parameters
	 */
	virtual void
	setRestartFilePath(const std::string& path) = 0;
};
} /* namespace options */
} /* namespace xolotl */
