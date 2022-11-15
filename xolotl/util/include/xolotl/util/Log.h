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
/**
 * @brief Serves as namespace for log-related functions as well as lazy
 * singleton for initializing logging
 */
class Log
{
public:
	/**
	 * Severity levels for log messages
	 */
	enum Level
	{
		debug = 0,
		extra,
		info,
		warning,
		error
	};

	/**
	 * Alias for severity logger type
	 */
	using LoggerType = boost::log::sources::severity_logger<Level>;

	Log();

    /**
     * @brief Set logging output threshold level
     */
    static void
    setLevelThreshold(Level level);

    /**
     * @brief Set logging output threshold level
     */
    static void
    setLevelThreshold(const std::string& levelString);

	/**
	 * @brief Get a severity logger
	 */
	static LoggerType&
	getLogger();

	/**
	 * @brief Flush the logging core
	 */
	static void
	flush();
};

/**
 * @relates Log
 * @brief Insert a Log::Level into a stream
 */
std::ostream&
operator<<(std::ostream& os, Log::Level level);

BOOST_LOG_INLINE_GLOBAL_LOGGER_DEFAULT(Logger, Log::LoggerType);

/**
 * @brief Customize a std::stringstream with a precision for floats
 */
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
