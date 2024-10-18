#pragma once

#include <xolotl/factory/Factory.h>

namespace xolotl
{
namespace interface
{
class IMaterialSubOptions;
} // namespace interface

namespace factory
{
namespace interface
{
class MaterialSubOptionsFactory :
	public Factory<MaterialSubOptionsFactory,
		::xolotl::interface::IMaterialSubOptions>
{
public:
	static std::string
	getFactoryName() noexcept
	{
		return "MaterialSubOptionsFactory";
	}

	static std::string
	getName(const options::IOptions& options)
	{
		return options.getMaterial();
	}
};
} // namespace interface
} // namespace factory
} // namespace xolotl
