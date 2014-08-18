package gov.ornl.xolotl.preprocessor;

import uk.co.flamingpenguin.jewel.cli.Option;

/**
 * This interface creates command line options for the parameters that are
 * needed to run Xolotl. If the options are not specified via the command line,
 * then the default values are used.
 */
public interface Arguments {

	/**
	 * This Option annotation corresponds to the '--maxHeSize' option which
	 * defines a default value of 8 and additionally provides a brief
	 * description of this option.
	 * 
	 * @param defaultValue
	 *            The default value for the maxium size of a Helium cluster in
	 *            the network if this option is not specified via the command
	 *            line
	 * @param description
	 *            Brief description of this option
	 */
	@Option(defaultValue = "8", description = "The maximum size of a Helium cluster in the network satisfying the condition 0 < maxHeSize < 9. (default = 8)")
	/**
	 * This operation produces the required command line option '--maxHeSize' which 
	 * takes a single integer value and is defined by the previous Option annotation
	 * 
	 * @return The maximum size of a Helium cluster in the network satisfying the condition 0 < maxHeSize < 9
	 */
	int getMaxHeSize();

	/**
	 * This Option annotation corresponds to the '--maxVSize' option which
	 * defines a default value of 29 and additionally provides a brief
	 * description of this option.
	 * 
	 * @param defaultValue
	 *            The default value for the maxium size of a vacancy cluster
	 *            in the network if this option is not specified via the
	 *            command line
	 * @param description
	 *            Brief description of this option
	 */
	@Option(defaultValue = "29", description = "The maximum size of a vacancy cluster in the network. (default = 29)")
	/**
	 * This operation produces the required command line option '--maxVSize' which 
	 * takes a single integer value and is defined by the previous Option annotation
	 * 
	 * @return The maxium size of a vacancy cluster in the network
	 */
	int getMaxVSize();

	/**
	 * This Option annotation corresponds to the '--maxISize' option which
	 * defines a default value of 6 and additionally provides a brief
	 * description of this option.
	 * 
	 * @param defaultValue
	 *            The default value for the maxium size of an interstitial
	 *            cluster in the network if this option is not specified via
	 *            the command line
	 * @param description
	 *            Brief description of this option
	 */
	@Option(defaultValue = "6", description = "The maximum size of an interstitial cluster in the network satisfying the condition 0 < maxISize < 7. (default = 6)")
	/**
	 * This operation produces the required command line option '--maxISize' which 
	 * takes a single integer value and is defined by the previous Option annotation
	 * 
	 * @return The maximum size of an interstitial cluster in the network satisfying the condition 0 < maxISize < 7
	 */
	int getMaxISize();

	/**
	 * This Option annotation corresponds to the '--startTemp' option which
	 * defines a default value of 1000 and additionally provides a brief
	 * description of this option.
	 * 
	 * @param defaultValue
	 *            The default, string, value for the constant starting
	 *            temperature (in Kelvin) if this option is not specified via
	 *            the command line
	 * @param description
	 *            Brief description of this option
	 */
	@Option(defaultValue = "1000", description = "The temperature (in Kelvin) will be the constant "
			+ "value specified (default = 1000)")
	/**
	 * This operation produces the required command line option '--startTemp' which 
	 * takes a single string value and is defined by the previous Option annotation
	 * 
	 * @return The constant starting temperature as a string
	 */
	String getStartTemp();

	/**
	 * This Option annotation corresponds to the '--perfHandler' option which
	 * defines the default value to use the standard performance handlers and
	 * additionally provides a brief description of the option.
	 * 
	 * @param defaultValue
	 *            The default performance handler that will be used if this
	 *            option is not specified via the command line
	 * @param description
	 *            Brief description of this option
	 */
	@Option(defaultValue = "std", description = "{dummy, std}  Which set of performance handlers to use (default = std)")
	/**
	 * This operation produces the required command line option '--perfHandler' which 
	 * takes a single string value and is defined by the previous Option annotation
	 * 
	 * @return The performance handler, as a string, that will be used
	 */
	String getPerfHandler();

	/**
	 * This Option annotation corresponds to the '--vizHandler' option which
	 * defines the default value to use the dummy visualization handlers and
	 * additionally provides a brief description of the option.
	 * 
	 * @param defaultValue
	 *            The default visualization handler that will be used if this
	 *            option is not specified via the command line
	 * @param description
	 *            Brief description of this option
	 */
	@Option(defaultValue = "dummy", description = "{dummy, std}  Which set of visualization handlers to use")
	/**
	 * This operation produces the required command line option '--vizHandler' which 
	 * takes a single string value and is defined by the previous Option annotation
	 * 
	 * @return The performance handler, as a string, that will be used
	 */
	String getVizHandler();

	/**
	 * This Option annotation corresponds to the '--petscArgs' option which
	 * defines the default value to be the single string of Petsc arguments and
	 * additionally provides a brief description of the option.
	 * 
	 * @param defaultValue
	 *            The single string of Petsc arguments that will be used if this
	 *            option is not specified via the command line
	 * @param description
	 *            Brief description of this option
	 */
	@Option(defaultValue = "-da_grid_x 10 -ts_final_time 1000 -ts_dt 1.0e-12 "
			+ "-ts_max_steps 3 -ts_adapt_dt_max 10 -ts_max_snes_failures 200 "
			+ "-pc_type fieldsplit -pc_fieldsplit_detect_coupling -fieldsplit_0_pc_type redundant "
			+ "-fieldsplit_1_pc_type sor -snes_monitor -ksp_monitor -ts_monitor", description = "List of arguments to be passed to PETSc")
	/**
	 * This operation produces the required command line option '--petscArgs' which 
	 * takes a single string value and is defined by the previous Option annotation
	 * 
	 * @return The single string of Petsc arguments
	 */
	String getPetscArgs();

	/**
	 * This Option annotation corresponds to the '--networkFile' option which
	 * defines a default network file name of networkInit.h5 and additionally
	 * provides a brief description of the option.
	 * 
	 * @param defaultValue
	 *            The default network file name that will be used if this option
	 *            is not specified via the command line
	 * @param description
	 *            Brief description of this option
	 */
	@Option(defaultValue = "networkInit.h5", description = "<fileName>  The HDF5 file containing the network (default = networkInit.h5)")
	/**
	 * This operation produces the required command line option '--networkFile' which 
	 * takes a single string value and is defined by the previous Option annotation
	 * 
	 * @return The name of the network file
	 */
	String getNetworkFile();

	/**
	 * This Option annotation corresponds to the optional '--material' option
	 * and provides a brief description of the option.
	 * 
	 * @param description
	 *            Brief description of this option
	 */
	@Option(description = "{W}  The optional option declaring which material will be used (Tungsten, W, is currently the only material Xolotl uses)")
	/**
	 * This operation produces the optional command line option '--material' which 
	 * takes a single string value and is defined by the previous Option annotation
	 * 
	 * @return The material
	 */
	String getMaterial();

	/**
	 * This operation makes the command line option '--material' optional.
	 * 
	 * @return Returns true if the option has been specified and false if it has
	 *         not
	 */
	boolean isMaterial();

	/**
	 * This Option annotation corresponds to the optional '--tempFile' option
	 * and provides a brief description of the option.
	 * 
	 * @param description
	 *            Brief description of this option
	 */
	@Option(description = "<tempFileName>  A temperature profile is given "
			+ "by the specified file, then linear interpolation is used to fit the data")
	/**
	 * This operation produces the optional command line option '--tempFile' which 
	 * takes a single string value and is defined by the previous Option annotation.
	 * NOTE:  This option should only be used when the user wishes to pass a file 
	 * containing a temperature profile to Xolotl.
	 * 
	 * @return The name of the temperature file
	 */
	String getTempFile();

	/**
	 * This operation makes the command line option '--tempFile' optional.
	 * 
	 * @return Returns true if the option has been specified and false if it has
	 *         not
	 */
	boolean isTempFile();

	/**
	 * This Option annotation corresponds to the optional '--heFlux' option and
	 * provides a brief description of the option.
	 * 
	 * @param description
	 *            Brief description of this option
	 */
	@Option(description = "This option allows the user to change the Helium flux by the factor specified (in nm).")
	/**
	 * This operation produces the optional command line option '--heFlux' which 
	 * takes a single string value and is defined by the previous Option annotation
	 * 
	 * @return The Helium flux
	 */
	String getHeFlux();

	/**
	 * This operation makes the command line option '--heFlux' optional.
	 * 
	 * @return Returns true if the option has been specified and false if it has
	 *         not
	 */
	boolean isHeFlux();

	/**
	 * This Option annotation corresponds to the optional '--maxHeFluence' option
	 * and provides a brief description of the option.
	 * 
	 * @param description
	 *            Brief description of this option
	 */
	@Option(description = "The maxium Helium fluence value (in nm) that will be used")
	// defaultValue = "10")
	/**
	 * This operation produces the optional command line option '--maxHeFluence' which 
	 * takes a single string value and is defined by the previous Option annotation
	 * 
	 * @return The maximum Helium fluence value
	 */
	String getMaxHeFluence();

	/**
	 * This operation makes the command line option '--maxHeFluence' optional.
	 * 
	 * @return Returns true if the option has been specified and false if it has
	 *         not
	 */
	boolean isMaxHeFluence();

	@Option(description = "<HDF5FileName>  Start the solver with the concentrations, time, "
			+ "and step size present in the HDF5 file")
	/**
	 * This operation produces the optional command line option '--checkpoint' which 
	 * takes a single string value and is defined by the previous Option annotation
	 * 
	 * @return The name of the file to read from to create the new HDF5 network file
	 */
	String getCheckpoint();

	/**
	 * This operation makes the command line option '--checkpoint' optional.
	 * 
	 * @return Returns true if the option has been specified and false if it has
	 *         not
	 */
	boolean isCheckpoint();

	/**
	 * This produces the command line arguments '--help' or '-h' either of which
	 * can be used to print usage help
	 * 
	 * @return True if the help message has been requested or false if it has
	 *         not
	 */
	@Option(helpRequest = true, description = "display help", shortName = "h")
	boolean getHelp();

}
