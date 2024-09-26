#pragma once

#include <boost/property_tree/ptree_fwd.hpp>

#include <xolotl/options/Options.h>

namespace xolotl
{
namespace options
{
namespace detail
{
class JSONElem;
class JSONElemVector;
} // namespace detail

class JSONOptions : public Options
{
public:
	JSONOptions();

	JSONOptions(const JSONOptions&);

	~JSONOptions();

	std::shared_ptr<IOptions>
	makeCopy() const override;

	void
	printHelp(std::ostream& = std::cout);

	void
	readParams(int argc, const char* argv[]) override;

private:
	detail::JSONElemVector
	defineHandlers();

	std::unique_ptr<boost::property_tree::iptree> _map;

	friend class detail::JSONElem;
};
} // namespace options
} // namespace xolotl
