#ifndef RADIUSMINSIZEOPTIONHANDLER_H
#define RADIUSMINSIZEOPTIONHANDLER_H

// Includes
#include "OptionHandler.h"
#include <algorithm>

namespace xolotlCore {

/**
 * RadiusMinSizeOptionHandler handles the number of dimensions option.
 */
class RadiusMinSizeOptionHandler: public OptionHandler {
public:

	/**
	 * The default constructor
	 */
	RadiusMinSizeOptionHandler() :
			OptionHandler("radiusSize",
					"radiusSize <minSize> ...              "
							"This option allows the user a minimum size for the computation for the average radius (default is 0).  \n") {
	}

	/**
	 * The destructor
	 */
	~RadiusMinSizeOptionHandler() {
	}

	/**
	 * This method will set the IOptions resoMinSize
	 * to the value given as the argument.
	 *
	 * @param opt The pointer to the option that will be modified.
	 * @param arg The minimum sizes.
	 */
	bool handler(IOptions *opt, const std::string& arg) {
		// Build an input stream from the argument string.
		xolotlCore::TokenizedLineReader<int> reader;
		auto argSS = std::make_shared<std::istringstream>(arg);
		reader.setInputStream(argSS);

		// Break the argument into tokens.
		auto tokens = reader.loadLine();

		// Create the array of sizes
		Array<int, 4> sizes;
		sizes.Init(0);

		// Set the values
		for (int i = 0; i < std::min((int) tokens.size(), 4); i++) {
			sizes[i] = tokens[i];
		}
		opt->setRadiusMinSizes(sizes);

		return true;
	}

};
//end class RadiusMinSizeOptionHandler

} /* namespace xolotlCore */

#endif
