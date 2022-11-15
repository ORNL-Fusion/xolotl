#include <algorithm>
#include <cctype>

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
#include <boost/phoenix.hpp>

#include <xolotl/util/Log.h>
#include <xolotl/util/MPIUtils.h>

namespace xolotl
{
namespace util
{
std::string
toUpper(const std::string& s)
{
	auto sret = s;
	std::transform(begin(sret), end(sret), begin(sret),
		[](unsigned char c) { return std::toupper(c); });
	return sret;
}

BOOST_LOG_ATTRIBUTE_KEYWORD(severity, "Severity", Log::Level)
BOOST_LOG_ATTRIBUTE_KEYWORD(process_id, "ProcessID",
	boost::log::attributes::current_process_id::value_type)

boost::log::attributes::current_process_id::value_type::native_type
getNativeProcessId(boost::log::value_ref<
	boost::log::attributes::current_process_id::value_type,
	tag::process_id> const& pid)
{
	if (pid) {
		return pid->native_id();
	}
	return 0;
}

using LevelType = std::underlying_type_t<Log::Level>;

static std::string levelStrings[] = {
	"DEBUG", "EXTRA", "INFO", "WARNING", "ERROR"};

const char*
toString(Log::Level level)
{
	auto l = static_cast<LevelType>(level);
	// if (l < sizeof(levelStrings) / sizeof(*levelStrings)) {
    if (l < sizeof(levelStrings)) {
		return levelStrings[level].c_str();
	}
	else {
		return nullptr;
	}
}

Log::Level
logLevelFromString(const std::string& levelStr)
{
	auto str = toUpper(levelStr);
	auto it = std::find(
		std::begin(levelStrings), std::end(levelStrings), str);
	if (it == std::end(levelStrings)) {
		throw std::runtime_error(
			"invalid logging output threshold level specified: " + levelStr);
	}
	return static_cast<Log::Level>(std::distance(std::begin(levelStrings), it));
}

std::ostream&
operator<<(std::ostream& os, Log::Level level)
{
	auto str = toString(level);
	if (str) {
		os << str;
	}
	else {
		os << static_cast<LevelType>(level);
	}

	return os;
}

Log::Log()
{
	namespace keywords = boost::log::keywords;
	namespace expr = boost::log::expressions;
	namespace attr = boost::log::attributes;

	auto core = boost::log::core::get();
#ifdef NDEBUG
	core->set_filter(severity > Log::debug);
#endif

	boost::log::add_common_attributes();
	core->add_global_attribute("Scope", attr::named_scope());

	// Print info and warning messages to stdout
	boost::log::add_console_log(
		std::cout, keywords::auto_flush = true, keywords::format = "%Message%")
		->set_filter(Log::info <= severity && severity < Log::error);
	// Print error messages to stderr
	boost::log::add_console_log(std::cerr, keywords::auto_flush = true,
		keywords::format = expr::stream << "[" << severity << "] "
										<< expr::message)
		->set_filter(severity >= Log::error);

	// Print all messages to log file with extra metadata
	bool mpiReady = mpiInitialized();
	std::string logFileBaseName = "xolotlOutput";
	std::string logFileExt = ".log";
	std::string mpiRankString = "";
	if (mpiReady) {
		auto rank = getMPIRank();
		auto rankStr = std::to_string(rank);
		logFileBaseName += "_" + rankStr;
		mpiRankString = "(R" + rankStr + ") ";
	}
	boost::log::add_file_log(keywords::file_name = logFileBaseName + logFileExt,
		keywords::format = expr::stream
			<< "("
			<< expr::format_date_time<boost::posix_time::ptime>(
				   "TimeStamp", "%Y-%m-%d %H:%M:%S")
			<< ") ="
			<< boost::phoenix::bind(&getNativeProcessId, process_id.or_none())
			<< "= " << mpiRankString << "[" << severity << "] "
			<< expr::format_named_scope("Scope", keywords::format = "{%n}")
			<< '\n'
			<< expr::message);

	if (!mpiReady) {
		BOOST_LOG_SEV(Logger::get(), Log::warning)
			<< "MPI rank will not be included in log messages because MPI has "
			   "not been initialized";
	}
}

void
Log::setLevelThreshold(Log::Level level)
{
	boost::log::core::get()->set_filter(severity >= level);
}

void
Log::setLevelThreshold(const std::string& levelString)
{
	setLevelThreshold(logLevelFromString(levelString));
}

Log::LoggerType&
Log::getLogger()
{
	static Log log;
	return Logger::get();
}

void
Log::flush()
{
	boost::log::core::get()->flush();
}
} // namespace util
} // namespace xolotl
