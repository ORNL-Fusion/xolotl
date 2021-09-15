#pragma once

#include <iomanip>
#include <sstream>

#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sources/record_ostream.hpp>
#include <boost/log/sources/severity_logger.hpp>

namespace xolotl
{
namespace util
{
class Log
{
public:
	enum Level
	{
		debug = 0,
		extra,
		info,
		warning,
		error
	};

	using LoggerType = boost::log::sources::severity_logger<Level>;

	Log();

	static LoggerType&
	getLogger();

	static void
	flush();
};

std::ostream&
operator<<(std::ostream& os, Log::Level level);

BOOST_LOG_INLINE_GLOBAL_LOGGER_DEFAULT(Logger, Log::LoggerType);

class StringStream : public std::stringstream
{
public:
	explicit StringStream(std::streamsize prec = 16)
	{
		this->precision(prec);
	}
};
} // namespace util
} // namespace xolotl

#define XOLOTL_LOG_SEV(level) \
	BOOST_LOG_SEV(::xolotl::util::Log::getLogger(), level) \
		<< std::setprecision(16)

#define XOLOTL_LOG XOLOTL_LOG_SEV(::xolotl::util::Log::info)

#define XOLOTL_LOG_DBG XOLOTL_LOG_SEV(::xolotl::util::Log::debug)
#define XOLOTL_LOG_XTRA XOLOTL_LOG_SEV(::xolotl::util::Log::extra)
#define XOLOTL_LOG_WARN XOLOTL_LOG_SEV(::xolotl::util::Log::warning)
#define XOLOTL_LOG_ERR XOLOTL_LOG_SEV(::xolotl::util::Log::error)
