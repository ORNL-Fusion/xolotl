#pragma once

namespace xolotl
{
namespace solver
{
namespace monitor
{
class IMonitor
{
public:
	virtual ~IMonitor()
	{
	}

	virtual void
	setup(int loop) = 0;

	virtual void
	keepFlux(std::vector<std::vector<std::vector<double>>>& nSurf,
		std::vector<std::vector<std::vector<double>>>& nBulk,
		std::vector<std::vector<std::vector<double>>>& surfFlux,
		std::vector<std::vector<std::vector<double>>>& bulkFlux) = 0;

	virtual void
	setFlux(std::vector<std::vector<std::vector<double>>>& nSurf,
		std::vector<std::vector<std::vector<double>>>& nBulk,
		std::vector<std::vector<std::vector<double>>>& surfFlux,
		std::vector<std::vector<std::vector<double>>>& bulkFlux) = 0;

	virtual void
	setExternalControlStep(std::size_t step) = 0;
};
} // namespace monitor
} // namespace solver
} // namespace xolotl
