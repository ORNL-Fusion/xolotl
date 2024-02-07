#pragma once

#include <xolotl/options/Options.h>

namespace xolotl
{
namespace options
{
class ConfOptions : public Options
{
public:
	ConfOptions();

	ConfOptions(const ConfOptions&);

	~ConfOptions();

	std::shared_ptr<IOptions>
	makeCopy() const override;

	void
	readParams(int argc, const char* argv[]) override;
};
} // namespace options
} // namespace xolotl
