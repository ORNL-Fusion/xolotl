#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>

#include <xolotl/util/Log.h>

namespace xolotl
{
namespace util
{
BOOST_LOG_GLOBAL_LOGGER_INIT(Logger, LoggerType)
{
	LoggerType lg;
	return lg;
}

void
initLogging()
{
#ifdef NDEBUG
	boost::log::core::get()->set_filter(
		boost::log::trivial::severity >= LogLevel::info);
#endif

	boost::log::add_common_attributes();

	namespace keywords = boost::log::keywords;
	namespace expr = boost::log::expressions;

	boost::log::add_console_log(std::cout, keywords::format = "%Message%")
		->set_filter(boost::log::trivial::severity < LogLevel::error);
	boost::log::add_console_log(
		std::cerr, keywords::format = "[%Severity%] %Message%")
		->set_filter(boost::log::trivial::severity >= LogLevel::error);

	boost::log::add_file_log(keywords::file_name = "xolotlOutput.log",
		keywords::format = "(%TimeStamp%)[%Severity%]: %Message%");
}

void
flushLogFile()
{
	boost::log::core::get()->flush();
}
} // namespace util
} // namespace xolotl
