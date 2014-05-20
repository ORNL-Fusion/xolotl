#ifndef LABELPROVIDER_H
#define LABELPROVIDER_H

// Includes
#include <Identifiable.h>
#include <string>

namespace xolotlViz {

/**
 * LabelProvider provides a series of labels to print on the plots.
 * One Plot must have only one LabelProvider.
 */
class LabelProvider: public xolotlCore::Identifiable {

private:

	/**
     * Declare the constructor as private to force the use of a name
	 */
    LabelProvider() : xolotlCore::Identifiable("unused") {};

public:

	/**
	 * The label of the X axis of the plot.
	 */
	std::string axis1Label = " ";

	/**
	 * The label of the Y axis of the plot.
	 */
	std::string axis2Label = " ";

	/**
	 * The label of the Z axis of the plot.
	 */
	std::string axis3Label = " ";

	/**
	 * The label for the time steps.
	 */
	std::string axis4Label = " ";

	/**
	 * Title label for the plot.
	 */
	std::string titleLabel = " ";

	/**
	 * Unit label for the plot.
	 */
	std::string unitLabel = " ";

	/**
	 * Time label for the plot.
	 */
	std::string timeLabel = " ";

	/**
	 * Time step label for the plot.
	 */
	std::string timeStepLabel = " ";

	/**
	 * The default constructor
	 */
	LabelProvider(std::string name) : xolotlCore::Identifiable(name) {};

	/**
	 * The destructor
	 */
	~LabelProvider() {};

};

//end class LabelProvider

} /* namespace xolotlViz */

#endif
