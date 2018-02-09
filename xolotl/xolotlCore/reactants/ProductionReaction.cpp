#include "ProductionReaction.h"

namespace xolotlCore {

std::ostream&
operator<<(std::ostream& os, const ProductionReaction& reaction) {

	os << '[' << "first: " << reaction.first << "; " << "second: "
			<< reaction.second << "; " << "rate: " << reaction.kConstant << ']';

	return os;
}

} // namespace xolotlCore

