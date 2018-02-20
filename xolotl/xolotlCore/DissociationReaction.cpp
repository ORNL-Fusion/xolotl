#include "DissociationReaction.h"

namespace xolotlCore {

std::ostream&
operator<<(std::ostream& os, const DissociationReaction& reaction) {

	os << '[' << "diss: " << reaction.dissociating << "; " << "first: "
			<< reaction.first << "; " << "second: " << reaction.second << "; "
			<< "rate: " << reaction.kConstant << "; " << "reverse: [";
	if (reaction.reverseReaction) {
		os << *(reaction.reverseReaction);
	} else {
		os << "none";
	}
	os << "]]";

	return os;
}

} // namespace xolotlCore

