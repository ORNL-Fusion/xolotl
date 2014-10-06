#ifndef OPTIONS_H
#define OPTIONS_H

#include "IOptions.h"
#include "optionhandlers/IOptionHandler.h"

namespace xolotlCore {

/**
 * Options realizes the IOptions interface.
 * All private members will be accessed through getters and setters.
 */
class Options : public IOptions {

protected:
	/**
	 * Map of options we support, keyed by option switch string
	 * (including leading dashes)
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
	 * The number of options that will be given to PETSc.
	 */
	int petscArgc;

	/**
	 * The pointer to the options that will be given to PETSc.
	 */
	char **petscArgv;

	/**
	 * The value of the step size for the spatial grid.
	 */
	double stepSize;

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
	 * Use the helium fluence option?
	 */
	bool heliumFluenceFlag;

	/**
	 * Value for the maximum fluence.
	 */
	double maxHeliumFluence;
	
	/**
	 * Use the helium flux option?
	 */
	bool heliumFluxFlag;

	/**
	 * Value for the  flux.
	 */
	double heliumFlux;

	/**
	 * Which type of performance infrastructure should we use?
	 */
	xolotlPerf::IHandlerRegistry::RegistryType perfRegistryType;

	/**
	 * Use the "standard" set of handlers for the visualization infrastructure?
	 */
	bool vizStandardHandlersFlag;

	/**
	 * Use the material option?
	 */
	bool materialFlag;

	/**
	 * Name of the material.
	 */
	std::string materialName;

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
    void readParams(int argc, char* argv[]);

    /**
     * Show our help message.
     * \see IOptions.h
     */
    void showHelp(std::ostream& os) const;

    /**
     * Should the program run after parsing the command line?
     * (Or did parsing the command line suggest the program
     * shouldn't/couldn't run?)
     * \see IOptions.h
     */
    bool shouldRun() const {return shouldRunFlag;}

    /**
     * Set the shouldRunFlag.
     * \see IOptions.h
     */
    void setShouldRunFlag(bool flag) {shouldRunFlag = flag;}

    /**
     * If program shouldn't run, what should its exit code be?
     * \see IOptions.h
     */
    int getExitCode() const {return exitCode;}

    /**
     * Set the value for the exit code
     * \see IOptions.h
     */
    void setExitCode(int code) {exitCode = code;}

    /**
     * Get the name of the network file.
     * \see IOptions.h
     */
    std::string getNetworkFilename() const {return networkFilename;}

    /**
     * Set the name of the network file.
     * \see IOptions.h
     */
    void setNetworkFilename(std::string name) {networkFilename = name;}

    /**
     * Get the Argc for PETSc.
     * \see IOptions.h
     */
    int getPetscArgc() const {return petscArgc;}

    /**
     * Set the Argc for PETSc.
     * \see IOptions.h
     */
    void setPetscArgc(int argc) {petscArgc = argc;}

    /**
     * Get the Argv for PETSc.
     * \see IOptions.h
     */
    char** getPetscArgv() const {return petscArgv;}

    /**
     * Set the Argv for PETSc.
     * \see IOptions.h
     */
    void setPetscArgv(char** argv) {petscArgv = argv;}

    /**
     * Get the value of the step size.
     * \see IOptions.h
     */
    double getStepSize() const {return stepSize;}

    /**
     * Set the value of the step size.
     * \see IOptions.h
     */
    void setStepSize(double value) {stepSize = value;}

    /**
     * Should we use const temperature handlers?
     * \see IOptions.h
     */
    bool useConstTemperatureHandlers() const {return constTempFlag;}

    /**
     * Set the constTempFlag.
     * \see IOptions.h
     */
    void setConstTempFlag(bool flag) {constTempFlag = flag;}

    /**
     * Obtain the value of the constant temperature to be used.
     * \see IOptions.h
     */
    double getConstTemperature() const {return constTemperature;}

    /**
     * Set the constant temperature.
     * \see IOptions.h
     */
    void setConstTemperature(double temp) {constTemperature = temp;}

    /**
     * Should we use temperature profile handlers?
     * \see IOptions.h
     */
    bool useTemperatureProfileHandlers() const {return tempProfileFlag;}

    /**
     * Set the tempProfileFlag.
     * \see IOptions.h
     */
    void setTempProfileFlag(bool flag) {tempProfileFlag = flag;}

    /**
     * Obtain the name of the file containing the temperature profile data.
     * \see IOptions.h
     */
    std::string getTempProfileFilename() const {return tempProfileFilename;}

    /**
     * Set the name of the profile file to use.
     * \see IOptions.h
     */
    void setTempProfileFilename(std::string name) {tempProfileFilename = name;}

    /**
     * Should we use the Helium fluence option?
	 * If false, it will not be used.
     * \see IOptions.h
     */
    bool useMaxHeliumFluence() const {return heliumFluenceFlag;}

    /**
     * Set the heliumFluenceFlag.
     * \see IOptions.h
     */
    void setHeliumFluenceFlag(bool flag) {heliumFluenceFlag = flag;}

    /**
     * Obtain the value of the Helium fluence to be used.
     * \see IOptions.h
     */
    double getMaxHeliumFluence() const {return maxHeliumFluence;}

    /**
     * Set the value for the maximum fluence to which we want to integrate.
     * \see IOptions.h
     */
    void setMaxHeliumFluence(double fluence) {maxHeliumFluence = fluence;}

    /**
     * Should we use the Helium flux option?
	 * If false, it will not be used.
     * \see IOptions.h
     */
    bool useHeliumFlux() const {return heliumFluxFlag;};

    /**
     * Set the heliumFluxFlag.
     * \see IOptions.h
     */
    void setHeliumFluxFlag(bool flag) {heliumFluxFlag = flag;}

    /**
     * Obtain the value of the Helium flux to be used.
     * \see IOptions.h
     */
    double getHeliumFlux() const {return heliumFlux;}

    /**
     * Set the value for the flux to use.
     * \see IOptions.h
     */
    void setHeliumFlux(double flux) {heliumFlux = flux;}

    /**
     * Which type of performance handlers should we use?
     * \see IOptions.h
     */
    xolotlPerf::IHandlerRegistry::RegistryType getPerfHandlerType(void) const { return perfRegistryType; }

    /**
     * Set the type of performance handlers to use.
     * \see IOptions.h
     */
    void setPerfHandlerType(xolotlPerf::IHandlerRegistry::RegistryType rtype) { perfRegistryType = rtype; }

    /**
     * Should we use the "standard" set of handlers for the visualization?
     * If false, use dummy (stub) handlers.
     * \see IOptions.h
     */
    bool useVizStandardHandlers() const {return vizStandardHandlersFlag;}

    /**
     * Set the vizStandardHandlersFlag.
     * \see IOptions.h
     */
    void setVizStandardHandlers(bool flag) {vizStandardHandlersFlag = flag;}

    /**
     * Should we use a specific material for the helium flux profile?
     * \see IOptions.h
     */
    bool useMaterial() const {return materialFlag;}

    /**
     * Set the materialFlag.
     * \see IOptions.h
     */
    void setMaterialFlag(bool flag) {materialFlag = flag;}

    /**
     * Obtain the name of the material to be used for the flux profile.
     * \see IOptions.h
     */
    std::string getMaterial() const {return materialName;}

    /**
     * Set the name of the material to be used for the flux profile.
     * \see IOptions.h
     */
    void setMaterial(std::string material) {materialName = material;}


};//end class Options

} /* namespace xolotlCore */

#endif // OPTIONS_H
