#ifndef IDENTIFIABLE_H
#define IDENTIFIABLE_H

#include <xolotl/core/IIdentifiable.h>

namespace xolotlCore {

/**
 * An Identifiable implements the IIdentifiable interface, so
 * that classes derived from Identifiable will be able to identify
 * themselves.
 */
class Identifiable: public virtual IIdentifiable {
private:

	std::string name;

public:

	Identifiable(const std::string& _name) :
			name(_name) {
	}

	virtual ~Identifiable(void) {
	}

	/**
	 * Obtain the object's given name.
	 */
	virtual std::string getName(void) const {
		return name;
	}
};

} // end namespace xolotlCore

#endif // IDENTIFIABLE_H
