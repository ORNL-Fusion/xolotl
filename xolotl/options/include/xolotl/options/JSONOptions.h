#pragma once

#include <boost/property_tree/ptree_fwd.hpp>

#include <xolotl/options/Options.h>

namespace xolotl
{
namespace options
{
class JSONOptions : public Options
{
public:
	JSONOptions();

	JSONOptions(const JSONOptions&);

	~JSONOptions();

	std::shared_ptr<IOptions>
	makeCopy() const override;

	void
	readParams(int argc, const char* argv[]) override;

private:
	std::unique_ptr<boost::property_tree::iptree> _map;
};
} // namespace options
} // namespace xolotl
