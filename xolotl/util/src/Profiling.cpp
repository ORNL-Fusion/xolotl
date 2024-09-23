#include <fstream>
#include <string>
#include <vector>

#include <xolotl/util/Profiling.h>

namespace xolotl
{
namespace util
{
static std::vector<std::string> regStack;

void
writeRegStack()
{
	auto fs = std::fstream("regions.txt", std::ios::out | std::ios::app);
	for (auto&& r : regStack) {
		fs << r << "    ";
	}
	fs << std::endl;
}

void
pushProfRegion(const std::string& label)
{
	Kokkos::Profiling::pushRegion(label);
	// regStack.push_back(label);
	// writeRegStack();
}

void
popProfRegion(const std::string& name)
{
	Kokkos::Profiling::popRegion();
	// if (!name.empty() && name != regStack.back()) {
	// 	auto fs = std::fstream("regions.txt", std::ios::out | std::ios::app);
	// 	fs << "Pop mismatch: {" << name << "}" << std::endl;
	// }
	// regStack.pop_back();
	// writeRegStack();
}
} // namespace util
} // namespace xolotl
