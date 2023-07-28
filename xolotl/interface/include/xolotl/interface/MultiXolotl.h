#pragma once

#include <memory>
#include <vector>

#include <xolotl/config.h>
#include <xolotl/interface/IXolotlInterface.h>
#include <xolotl/options/IOptions.h>

namespace xolotl
{
namespace interface
{
class ComputeContext;
class XolotlInterface;

class MultiXolotl : public IXolotlInterface
{
public:
	MultiXolotl(const std::shared_ptr<ComputeContext>& context,
		const std::shared_ptr<options::IOptions>& opts);

	virtual ~MultiXolotl();

	void
	solveXolotl() override;

private:
	std::shared_ptr<ComputeContext> _computeContext;
	std::shared_ptr<options::IOptions> _options;
	std::unique_ptr<XolotlInterface> _primaryInstance;
	std::vector<std::unique_ptr<XolotlInterface>> _subInstances;
	std::vector<IdType> _subDOFs;
	double _currentTime{};
	double _maxDt{};
	IdType _localXM{1};
};

// class XolotlNetworkProblem
// {
// public:
// 	XolotlNetworkProblem(const InputParameters& params);

// 	~XolotlNetworkProblem();

// 	void
// 	externalSolve();
// 	bool
// 	converged();

// 	// Methods for restart
// 	void
// 	saveState();
// 	void
// 	setState();

// private:
// 	/// The path to the input file for Xolotl
// 	FileName _network_xolotl_filename;
// 	std::vector<FileName> _subnetwork_xolotl_filenames;
// 	std::shared_ptr<XolotlInterface> _networkInterface;
// 	std::vector<std::shared_ptr<XolotlInterface>> _subInterfaces;
// 	std::vector<xolotl::IdType> _subDOFs;
// 	Real& _current_time;
// 	Real _max_dt;
// 	xolotl::IdType _localXS;
// 	xolotl::IdType _localXM;

// 	// Variables for restart
// 	Real& _current_dt;
// 	Real& _previous_time;
// 	std::vector<
// 		std::vector<std::vector<std::vector<std::pair<xolotl::IdType, Real>>>>>&
// 		_conc_vector;
// };
} // namespace interface
} // namespace xolotl
