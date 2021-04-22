#pragma once

#include <array>
#include <string>

namespace xolotl
{
namespace test
{
template <int NArgs>
struct CommandLine
{
	CommandLine() = default;

	CommandLine(const std::array<std::string, NArgs>& args) : _args(args)
	{
		for (int i = 0; i < argc; ++i) {
			argv[i] = _args[i].c_str();
		}
	}

	int argc{NArgs};
	const char* argv[NArgs + 1]{nullptr};

private:
	const std::array<std::string, NArgs> _args;
};
} // namespace test
} // namespace xolotl
