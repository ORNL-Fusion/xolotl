#pragma once

#include <memory>
#include <vector>

#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace interface
{
class IMaterialSubOptions
{
public:
	virtual ~IMaterialSubOptions()
	{
	}

	virtual std::vector<std::shared_ptr<options::IOptions>>
	getSubOptions() const = 0;
};
} // namespace interface
} // namespace xolotl
