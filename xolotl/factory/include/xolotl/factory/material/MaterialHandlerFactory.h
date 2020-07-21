#pragma once

#include <xolotl/factory/Factory.h>

namespace xolotl
{
namespace core
{
namespace material
{
class IMaterialHandler;
}
} // namespace core

namespace factory
{
namespace material
{
class MaterialHandlerFactory :
    public Factory<MaterialHandlerFactory, core::material::IMaterialHandler>
{
public:
    static std::string
    getFactoryName() noexcept
    {
        return "MaterialHandlerFactory";
    }

    static std::string
    getName(const options::Options& options)
    {
        return options.getMaterial();
    }
};
} // namespace material
} // namespace factory
} // namespace xolotl
