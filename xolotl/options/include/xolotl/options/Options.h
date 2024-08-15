#pragma once

#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace options
{
/**
 * Options realizes the IOptions interface.
 * All private members will be accessed through getters and setters.
 */
class Options : public IOptions
{
protected:
	/**
	 * The name of the file to use for checkpoint
	 */
	std::string checkpointFile;

	/**
	 * The name of the file to use for restart
	 */
	std::string restartFile;

	/**
	 * The options that will be given to PETSc.
	 */
	std::string petscArg;

	/**
	 * Name of the temperature handler to use
	 */
	std::string tempHandlerName;

	/**
	 * Temperature parameter values
	 */
	std::array<double, 2> tempParam;

	/**
	 * Value for the temperature grid power.
	 */
	double tempGridPower;

	/**
	 * Name of the input temperature profile file.
	 */
	std::string tempProfileFilename;

	/**
	 * Use the flux amplitude option?
	 */
	bool fluxFlag;

	/**
	 * Value for the flux.
	 */
	double fluxAmplitude;

	/**
	 * Use a time profile for the flux?
	 */
	bool fluxTimeProfileFlag;

	/**
	 * Name of the input time profile file for the flux.
	 */
	std::string fluxTimeProfileFilePath;

	/**
	 * Value for the helium portion.
	 */
	double hePPM;

	/**
	 * Name of the perf handler
	 */
	std::string perfHandlerName;

	/**
	 * Output performance report to YAML file?
	 */
	bool perfOutputYAMLFlag;

	/**
	 * Enable multiple xolotl instances
	 */
	bool subnetworksFlag;

	/**
	 * Initial coupling timestep
	 */
	double initialTimeStep;

	/**
	 * Maximum coupling timestep
	 */
	double maxTimeStep;

	/**
	 * Time step growth factor
	 */
	double timeStepGrowthFactor;

	/**
	 * Coupling start time
	 */
	double startTime;

	/**
	 * Coupling end time
	 */
	double endTime;

	/**
	 * Maximum coupling time steps
	 */
	IdType numTimeSteps;

	/**
	 * Name of the viz handler
	 */
	std::string vizHandlerName;

	/**
	 * Name of the material.
	 */
	std::string materialName;

	/**
	 * Value for initial concentrations.
	 */
	std::string initialConcentration;

	/**
	 * Value of the electronic stopping power.
	 */
	double zeta;

	/**
	 * The location of the interface between two materials.
	 */
	double interfaceLocation;

	/**
	 * Number of dimensions for the simulation.
	 */
	int dimensionNumber;

	/**
	 * Name of the grid type to use
	 */
	std::string gridTypeName;

	/**
	 * Grid parameter values
	 */
	std::array<double, 6> gridParam;

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
	 * Network parameters
	 */
	std::vector<IdType> networkParams;

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
	 * Maximum number of pure V
	 */
	int maxPureV;

	/**
	 * Maximum number of I
	 */
	int maxI;

	/**
	 * The boundary condition on the given side of the grid.
	 */
	int leftBoundary, rightBoundary, bottomBoundary, topBoundary, frontBoundary,
		backBoundary;

	/**
	 * String of the list of wanted BC in X.
	 */
	std::string xBC;

	/**
	 * Portion of heat lost in the bulk.
	 */
	double heatLossPortion;

	/**
	 * Depth for the bubble bursting in nm.
	 */
	double burstingDepth;

	/**
	 * Factor used in bursting probability.
	 */
	double burstingFactor;

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
	 * Average radius computation minimum size
	 */
	std::vector<size_t> radiusMinSizes;

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
	 * Reflect the fact that interstitial clusters have a larger surrounding
	 * strain field.
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
	 * HeV ration, how many He per V are allowed
	 */
	double heVRatio;

	/**
	 * Migration energy above which the diffusion will be ignored
	 */
	double migrationThreshold;

	/**
	 * The path to the custom flux profile file
	 */
	fs::path fluxDepthProfileFilePath;

	/**
	 * The path to the reaction rates file
	 */
	fs::path reactionFilePath;

	/**
	 * Value of the basal portion.
	 */
	double basalPortion;

	/**
	 * Transition size, for instance from pyramic to c-loops.
	 */
	int transitionSize;

	/**
	 * Value of the cascade dose.
	 */
	double cascadeDose;

	/**
	 * Value of the remaining cascade efficiency.
	 */
	double cascadeEfficiency;

public:
	/**
	 * The constructor.
	 */
	Options();

	/**
	 * The destructor.
	 */
	~Options();

	void
	printAll(std::ostream& os) const override;

	/**
	 * \see IOptions.h
	 */
	std::string
	getCheckpointFilePath() const override
	{
		return checkpointFile;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getRestartFilePath() const override
	{
		return restartFile;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getSolverName() const override
	{
		return "PETSc";
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getPetscArg() const override
	{
		return petscArg;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getTempHandlerName() const override
	{
		return tempHandlerName;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getTempParam(std::size_t i = 0) const override
	{
		return tempParam[i];
	}

	/**
	 * \see IOptions.h
	 */
	double
	getTempGridPower() const override
	{
		return tempGridPower;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getTempProfileFilename() const override
	{
		return tempProfileFilename;
	}

	/**
	 * \see IOptions.h
	 */
	bool
	useFluxAmplitude() const override
	{
		return fluxFlag;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getFluxAmplitude() const override
	{
		return fluxAmplitude;
	}

	/**
	 * \see IOptions.h
	 */
	bool
	useFluxTimeProfile() const override
	{
		return fluxTimeProfileFlag;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getFluxTimeProfileFilePath() const override
	{
		return fluxTimeProfileFilePath;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getHePPM() const override
	{
		return hePPM;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getPerfHandlerName() const override
	{
		return perfHandlerName;
	}

	/**
	 * \see IOptions.h
	 */
	bool
	usePerfOutputYAML() const override
	{
		return perfOutputYAMLFlag;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getVizHandlerName() const override
	{
		return vizHandlerName;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getMaterial() const override
	{
		return materialName;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getInitialConcentration() const override
	{
		return initialConcentration;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getZeta() const override
	{
		return zeta;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getDimensionNumber() const override
	{
		return dimensionNumber;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getInterfaceLocation() const override
	{
		return interfaceLocation;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getGridTypeName() const override
	{
		return gridTypeName;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getGridParam(std::size_t i = 0) const override
	{
		return gridParam[i];
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getGridFilename() const override
	{
		return gridFilename;
	}

	/**
	 * \see IOptions.h
	 */
	const std::map<std::string, bool>&
	getProcesses() const override
	{
		return processMap;
	}

	void
	addProcess(const std::string& processKey) override;

	/**
	 * \see IOptions.h
	 */
	bool
	useSubnetworks() const override
	{
		return subnetworksFlag;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getInitialTimeStep() const override
	{
		return initialTimeStep;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getMaxTimeStep() const override
	{
		return maxTimeStep;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getTimeStepGrowthFactor() const override
	{
		return timeStepGrowthFactor;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getStartTime() const override
	{
		return startTime;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getEndTime() const override
	{
		return endTime;
	}

	/**
	 * \see IOptions.h
	 */
	virtual IdType
	getNumberOfTimeSteps() const override
	{
		return numTimeSteps;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getGbString() const override
	{
		return gbList;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getGroupingMin() const override
	{
		return groupingMin;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getGroupingWidthA() const override
	{
		return groupingWidthA;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getGroupingWidthB() const override
	{
		return groupingWidthB;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getSputteringYield() const override
	{
		return sputteringYield;
	}

	/**
	 * \see IOptions.h
	 */
	bool
	useHDF5() const override
	{
		return useHDF5Flag;
	}

	/**
	 * \see IOptions.h
	 */
	const std::vector<IdType>&
	getNetworkParameters() const override
	{
		return networkParams;
	}

	/**
	 * \see IOptions.h
	 */
	void
	setNetworkParameters(const std::vector<IdType>& params) override;

	/**
	 * \see IOptions.h
	 */
	int
	getMaxImpurity() const override
	{
		return maxImpurity;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getMaxD() const override
	{
		return maxD;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getMaxT() const override
	{
		return maxT;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getMaxV() const override
	{
		return maxV;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getMaxPureV() const override
	{
		return maxPureV;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getMaxI() const override
	{
		return maxI;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getLeftBoundary() const override
	{
		return leftBoundary;
	}
	int
	getRightBoundary() const override
	{
		return rightBoundary;
	}
	int
	getBottomBoundary() const override
	{
		return bottomBoundary;
	}
	int
	getTopBoundary() const override
	{
		return topBoundary;
	}
	int
	getFrontBoundary() const override
	{
		return frontBoundary;
	}
	int
	getBackBoundary() const override
	{
		return backBoundary;
	}

	/**
	 * \see IOptions.h
	 */
	std::string
	getBCString() const override
	{
		return xBC;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getHeatLossPortion() const override
	{
		return heatLossPortion;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getBurstingDepth() const override
	{
		return burstingDepth;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getBurstingFactor() const override
	{
		return burstingFactor;
	}

	/**
	 * \see IOptions.h
	 */
	void
	setRNGSeed(unsigned int s) override
	{
		rngSeed = s;
		rngUseSeed = true;
	}

	/**
	 * \see IOptions.h
	 */
	std::tuple<bool, unsigned int>
	getRNGSeed(void) const override
	{
		return std::make_tuple(rngUseSeed, rngSeed);
	}

	/**
	 * \see IOptions.h
	 */
	bool
	printRNGSeed(void) const override
	{
		return rngPrintSeed;
	}

	/**
	 * \see IOptions.h
	 */
	virtual std::vector<size_t>
	getRadiusMinSizes() const override
	{
		return radiusMinSizes;
	}

	/**
	 * \see IOptions.h
	 */
	double
	getDensity() const override
	{
		return density;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getPulseTime() const override
	{
		return pulseTime;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getPulseProportion() const override
	{
		return pulseProportion;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getLatticeParameter() const override
	{
		return latticeParameter;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getImpurityRadius() const override
	{
		return impurityRadius;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getBiasFactor() const override
	{
		return biasFactor;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getHydrogenFactor() const override
	{
		return hydrogenFactor;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getXenonDiffusivity() const override
	{
		return xenonDiffusivity;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getFissionYield() const override
	{
		return fissionYield;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getHeVRatio() const override
	{
		return heVRatio;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getMigrationThreshold() const override
	{
		return migrationThreshold;
	}

	/**
	 * \see IOptions.h
	 */
	virtual std::string
	getFluxDepthProfileFilePath() const override
	{
		return fluxDepthProfileFilePath.string();
	}

	/**
	 * \see IOptions.h
	 */
	virtual std::string
	getReactionFilePath() const override
	{
		return reactionFilePath.string();
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getBasalPortion() const override
	{
		return basalPortion;
	}

	/**
	 * \see IOptions.h
	 */
	int
	getTransitionSize() const override
	{
		return transitionSize;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getCascadeDose() const override
	{
		return cascadeDose;
	}

	/**
	 * \see IOptions.h
	 */
	virtual double
	getCascadeEfficiency() const override
	{
		return cascadeEfficiency;
	}

protected:
	/**
	 * \see IOptions.h
	 */
	void
	setCheckpointFilePath(const std::string& path) override
	{
		checkpointFile = path;
	}

	/**
	 * \see IOptions.h
	 */
	void
	setRestartFilePath(const std::string& path) override
	{
		restartFile = path;
	}

	void
	setNetworkParameters(const std::string& paramStr);

	void
	setTempParam(const std::vector<double>& params);

	void
	setTempParam(const std::string& paramStr);

	void
	checkTempProfileFilename() const;

	void
	checkFluxTimeProfileFilename() const;

	void
	checkPerfHandlerName() const;

	void
	checkVizHandlerName() const;

	void
	setGridParam(const std::vector<double>& params);

	void
	setGridParam(const std::string& paramStr);

	void
	checkGridFilename() const;

	void
	setRadiusMinSizes(const std::vector<int>& params);

	void
	setRadiusMinSizes(const std::string& paramStr);

	void
	setBoundaries(const std::vector<int>& params);

	void
	setBoundaries(const std::string& paramStr);

	void
	processRNGParam(const std::string& paramStr);

	void
	setProcesses(const std::vector<std::string>& processList);

	void
	setProcesses(const std::string& processList);

	void
	setCouplingTimeStepParams(const std::vector<double>& params);

	void
	setCouplingTimeStepParams(const std::string& paramString);

	void
	setPulseParams(const std::string& paramStr);

	void
	setGroupingParams(const std::vector<int>& params);

	void
	setGroupingParams(const std::string& paramString);

	void
	appendPetscArg(const std::string& arg);
};

std::shared_ptr<IOptions>
createOptions(int argc, const char* argv[]);
} /* namespace options */
} /* namespace xolotl */
