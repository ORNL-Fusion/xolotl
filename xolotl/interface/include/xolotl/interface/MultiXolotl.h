#pragma once

#include <memory>

#include <xolotl/interface/IXolotlInterface.h>

namespace xolotl
{
namespace options
{
class IOptions;
}
namespace interface
{
class ComputeContext;

class MultiXolotl : public IXolotlInterface
{
public:
	MultiXolotl(const std::shared_ptr<ComputeContext>& context,
		const std::shared_ptr<options::IOptions>& opts);

	virtual ~MultiXolotl();

	void
	solveXolotl() override;
};
} // namespace interface
} // namespace xolotl
