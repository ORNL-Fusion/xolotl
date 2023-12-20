#pragma once

#include <string>

#include <xolotl/core/temperature/TemperatureHandler.h>
#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace core
{
namespace temperature
{
/**
 * This class realizes the TemperatureHandler, it is responsible for the
 * handling of a temperature changing with time.
 */
class ProfileHandler : public TemperatureHandler
{
private:
	/**
	 * The name of the file were the profile is stored.
	 */
	std::string tempFile;

	/**
	 * Vector to hold the time read from the input
	 * temperature file.
	 */
	std::vector<double> time;

	/**
	 * Vector to hold the temperature read from the input
	 * temperature file.
	 */
	std::vector<double> temp;

public:
	ProfileHandler() = delete;

	/**
	 * The constructor.
	 *
	 * @param profileFileName The name of the profile file
	 */
	ProfileHandler(const std::string& profileFileName);

	/**
	 * Construct from options
	 */
	ProfileHandler(const options::IOptions& options);

	/**
	 * The destructor.
	 */
	virtual ~ProfileHandler();

	/**
	 * \see ITemperatureHandler.h
	 */
	void
	initialize(const int dof) override;

	/**
	 * This operation returns the temperature at the given position
	 * and time.
	 * It linearly interpolates the data read from the input
	 * temperature file.
	 *
	 * \see ITemperatureHandler.h
	 */
	double
	getTemperature(
		const plsm::SpaceVector<double, 3>&, double currentTime) const override;
};
// end class ProfileHandler

} // namespace temperature
} // namespace core
} // namespace xolotl
