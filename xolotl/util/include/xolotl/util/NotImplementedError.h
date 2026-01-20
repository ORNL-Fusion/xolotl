#include <stdexcept>

namespace xolotl
{
namespace util
{
class NotImplementedError : public std::runtime_error
{
public:
	NotImplementedError() : std::runtime_error("This has not been implemented.")
	{
	}

	using std::runtime_error::runtime_error;
};
} // namespace util
} // namespace xolotl
