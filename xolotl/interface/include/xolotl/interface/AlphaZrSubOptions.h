#pragma once

#include <xolotl/interface/IMaterialSubOptions.h>

namespace xolotl
{
namespace interface
{
class AlphaZrSubOptions : public IMaterialSubOptions
{
public:
	AlphaZrSubOptions(const options::IOptions& options);

	std::vector<std::shared_ptr<options::IOptions>>
	getSubOptions() const override;

private:
	const options::IOptions& _options;
};
} // namespace interface
} // namespace xolotl
