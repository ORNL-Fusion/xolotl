#include <xolotl/util/StreamUtils.h>

namespace xolotl
{
namespace util
{
std::stringstream
stripComments(std::istream& is)
{
	std::stringstream ss;
	std::string line;
	bool cCommenting = false;
	while (std::getline(is, line)) {
		for (auto it = begin(line); it != end(line); ++it) {
			if (cCommenting) {
				if (*it == '*' && next(it) != end(line) && *next(it) == '/') {
					cCommenting = false;
					++it;
				}
				continue;
			}
			if (*it == '#') {
				break;
			}
			if (*it == '/') {
				if (next(it) != end(line) && *next(it) == '/') {
					break;
				}
				if (next(it) != end(line) && *next(it) == '*') {
					cCommenting = true;
					++it;
					continue;
				}
			}
			ss << *it;
		}
	}
	return ss;
}
}
}
