#pragma once

#include <iomanip>
#include <sstream>

#include <boost/log/attributes/named_scope.hpp>
#include <boost/log/sources/global_logger_storage.hpp>
#include <boost/log/sources/severity_logger.hpp>
#include <boost/log/sources/record_ostream.hpp>

namespace xolotl
{
namespace util
{
enum LogLevel
{
	debug,
	extra,
	info,
	warning,
	error
};

const char*
toString(LogLevel level);

std::ostream&
operator<<(std::ostream& os, LogLevel level);

using LoggerType = boost::log::sources::severity_logger<LogLevel>;

BOOST_LOG_GLOBAL_LOGGER(Logger, LoggerType);

inline decltype(auto)
getLogger()
{
	return Logger::get();
}

void
initLogging();

void
flushLogFile();

std::string
getLogFileName();

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
	BOOST_LOG_SEV(::xolotl::util::getLogger(), level) << std::setprecision(16)

#define XOLOTL_LOG XOLOTL_LOG_SEV(::xolotl::util::LogLevel::info)

#define XOLOTL_LOG_DBG XOLOTL_LOG_SEV(::xolotl::util::LogLevel::debug)
#define XOLOTL_LOG_XTRA XOLOTL_LOG_SEV(::xolotl::util::LogLevel::extra)
#define XOLOTL_LOG_WARN XOLOTL_LOG_SEV(::xolotl::util::LogLevel::warning)
#define XOLOTL_LOG_ERR XOLOTL_LOG_SEV(::xolotl::util::LogLevel::error)
