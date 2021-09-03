#include <boost/log/attributes.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>
#include <boost/log/expressions/formatters/named_scope.hpp>
#include <boost/log/sinks/text_file_backend.hpp>
#include <boost/log/support/date_time.hpp>
#include <boost/log/support/exception.hpp>
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>

#include <xolotl/util/Log.h>

namespace xolotl
{
namespace util
{
BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", LogLevel)

const char*
toString(LogLevel level)
{
	static const char* strings[] = {
		"DEBUG", "EXTRA", "INFO", "WARNING", "ERROR"};

	auto l = static_cast<std::size_t>(level);
	if (l < sizeof(strings) / sizeof(*strings)) {
		return strings[level];
	}
	else {
		return nullptr;
	}
}

std::ostream&
operator<<(std::ostream& os, LogLevel level)
{
	auto str = toString(level);
	if (str) {
		os << str;
	}
	else {
		os << static_cast<int>(level);
	}

	return os;
}

BOOST_LOG_GLOBAL_LOGGER_INIT(Logger, LoggerType)
{
	LoggerType lg;
	return lg;
}

void
initLogging()
{
	namespace keywords = boost::log::keywords;
	namespace expr = boost::log::expressions;
	namespace attr = boost::log::attributes;

	auto core = boost::log::core::get();
#ifdef NDEBUG
	core->set_filter(severity > LogLevel::debug);
#endif

	boost::log::add_common_attributes();
	core->add_global_attribute("Scope", attr::named_scope());

	boost::log::add_console_log(std::cout, keywords::format = "%Message%")
		->set_filter(LogLevel::info <= severity && severity < LogLevel::error);
	boost::log::add_console_log(std::cerr,
		keywords::format = expr::stream << "[" << severity << "] "
										<< expr::message)
		->set_filter(severity >= LogLevel::error);

	boost::log::add_file_log(keywords::file_name = "xolotlOutput.log",
		keywords::format = expr::stream
			<< "("
			<< expr::format_date_time<boost::posix_time::ptime>(
				   "TimeStamp", "%Y-%m-%d %H:%M:%S")
			<< ")[" << severity << "] " << expr::message
			<< expr::format_named_scope("Scope", keywords::format = "{%n}"));
}

void
flushLogFile()
{
	boost::log::core::get()->flush();
}
} // namespace util
} // namespace xolotl
