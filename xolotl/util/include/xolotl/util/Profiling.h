#pragma once

#include <Kokkos_Core.hpp>

namespace xolotl
{
namespace util
{
// inline
void
pushProfRegion(const std::string& label);
// {
// 	Kokkos::Profiling::pushRegion(label);
// }

// inline
void
popProfRegion(const std::string& name = "");
// {
// 	Kokkos::Profiling::popRegion();
// }

struct ProfileRegion
{
	ProfileRegion() = default;

	ProfileRegion(const std::string& label)
	{
		// pushProfRegion(label);
		Kokkos::Profiling::pushRegion(label);
	}

	~ProfileRegion()
	{
		if (!active) {
			// popProfRegion();
			Kokkos::Profiling::popRegion();
		}
	}

	bool active{false};
};
} // namespace util
} // namespace xolotl

#define XOLOTL_PROF_REGION(name) \
	auto __region = ::xolotl::util::ProfileRegion(name);
