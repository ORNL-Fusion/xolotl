#include <cassert>
#include <iomanip>
#include <sstream>
#include <string_view>

#include <xolotl/options/detail/JSONElem.h>
#include <xolotl/util/Tokenizer.h>

namespace xolotl::options::detail
{

inline const std::string&
typeString(const JSONElem& elem)
{
	static const std::string labelArray[] = {
		"BOOL", "INT", "STR", "REAL", "LIST{INT}", "LIST{STR}", "LIST{REAL}"};
	return labelArray[static_cast<int>(elem.type())];
}

std::stringstream
formatFirstColumn(const JSONElem& elem)
{
	std::stringstream ss;
	ss << "  " << elem.name();
	if (elem.type() != JSONElem::Type::invalid) {
		ss << " (" << typeString(elem) << ")";
	}
	return ss;
}

namespace
{
// This is stolen (and simplified) from Boost.ProgramOptions
inline constexpr unsigned defaultLineLength = 80;

unsigned
getLineLength()
{
	return defaultLineLength;
}

unsigned
getMinHelpLength()
{
	return defaultLineLength / 2;
}

unsigned
getElemColumnWidth(const std::deque<JSONElem>& elems, const unsigned lineLength,
	const unsigned minHelpLength)
{
	const unsigned startOfHelpColumn = lineLength - minHelpLength;

	// Find the maximum width of the option column
	unsigned width{23};
	for (auto&& elem : elems) {
		auto ss = formatFirstColumn(elem);
		width = std::max(width, static_cast<unsigned>(ss.str().size()));
	}

	// If first column is longer than the start of the help column,
	// we'll go to a new line
	width = std::min(width, startOfHelpColumn - 1);

	// add an additional space to improve readability
	++width;
	return width;
}

void
leftPad(std::ostream& os, unsigned width)
{
	os << std::setfill(' ') << std::setw(width) << "";
}

// NOTE: a constructor of this type is provided in C++20
template <typename It, typename End>
auto
stringView(It start, End stop)
{
	return std::string_view(&(*start), std::distance(start, stop));
}

void
printParagraph(std::ostream& os, std::string_view par, unsigned indent,
	unsigned lineLength)
{
	// Through reminder of this function, 'lineLength' will
	// be the length available for characters, not including
	// indent.
	assert(indent < lineLength);
	lineLength -= indent;

	if (par.size() < lineLength) {
		os << par;
	}
	else {
		const auto parEnd = par.cend();
		for (auto lineBegin = par.cbegin(); lineBegin != parEnd;) {
			// If line starts with space, but second character
			// is not space, remove the leading space.
			// We don't remove double spaces because those
			// might be intentianal.
			if ((*lineBegin == ' ') &&
				((lineBegin + 1 < parEnd) && (*(lineBegin + 1) != ' '))) {
				++lineBegin;
			}

			unsigned remaining =
				static_cast<unsigned>(std::distance(lineBegin, parEnd));
			auto lineEnd =
				lineBegin + ((remaining < lineLength) ? remaining : lineLength);

			// prevent chopped words
			// Is lineEnd between two non-space characters?
			if ((*(lineEnd - 1) != ' ') &&
				((lineEnd < parEnd) && (*lineEnd != ' '))) {
				// find last ' ' in the second half of the current paragraph
				// line
				auto line = stringView(lineBegin, lineEnd);
				auto lastSpace = line.find_last_of(' ');

				if (lastSpace != std::string_view::npos) {
					// is lastSpace within the second half of the current line
					if ((line.size() - lastSpace) < (lineLength / 2)) {
						lineEnd = lineBegin + lastSpace;
					}
				}
			}

			// write line to stream
			os << stringView(lineBegin, lineEnd);

			// more lines to follow?
			if (lineEnd != parEnd) {
				os << '\n';
				leftPad(os, indent);
			}

			lineBegin = lineEnd;
		}
	}
}

void
printElemHelp(std::ostream& os, const std::string& desc, unsigned firstColWidth,
	unsigned lineLength)
{
	// we need to use one char less per line to work correctly if actual
	// console has longer lines
	assert(lineLength > 1);
	if (lineLength > 1) {
		--lineLength;
	}

	// lineLength must be larger than firstColWidth
	// this assert may fail due to user error or environment conditions!
	assert(lineLength > firstColWidth);

	auto paragraphs = util::Tokenizer<std::string>{desc, "\n"}();

	unsigned padSize = 0;
	for (auto&& par : paragraphs) {
		leftPad(os, padSize);
		printParagraph(os, par, firstColWidth, lineLength);
		os << '\n';
		padSize = firstColWidth;
	}
}

void
printElem(std::ostream& os, const JSONElem& elem, unsigned firstColWidth,
	unsigned lineLength)
{
	auto ss = formatFirstColumn(elem);
	os << ss.str();

	if (!elem.help().empty()) {
		unsigned padSize =
			firstColWidth - static_cast<unsigned>(ss.str().size());
		if (ss.str().size() >= firstColWidth) {
			os.put('\n'); // first column is too long, new line for help
			padSize = firstColWidth;
		}
		leftPad(os, padSize);

		printElemHelp(os, elem.help(), firstColWidth, lineLength);
	}
}
} // namespace

void
JSONElemVector::print(std::ostream& os) const
{
	auto lineLength = getLineLength();
	auto minHelpLength = getMinHelpLength();
	auto width = getElemColumnWidth(_elems, lineLength, minHelpLength);

	for (auto&& elem : _elems) {
		printElem(os, elem, width, lineLength);
		os << "\n";
	}
}

void
JSONElemVector::processParams()
{
	for (auto&& elem : _elems) {
		elem();
	}
}
} // namespace xolotl::options::detail
