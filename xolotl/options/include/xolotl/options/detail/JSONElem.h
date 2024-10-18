#pragma once

#include <deque>
#include <functional>
#include <string>
#include <utility>

namespace xolotl
{
namespace options
{
namespace detail
{
class JSONElem
{
public:
	enum class Type
	{
		boolean,
		integer,
		string,
		real,
		list_integer,
		list_string,
		list_real,
		invalid
	};

	JSONElem() = delete;

	template <typename F>
	JSONElem(const std::string& name, const std::string& help, F func) :
		_name(name),
		_help(help),
		_func(func)
	{
	}

	template <typename F>
	JSONElem(
		const std::string& name, Type type, const std::string& help, F func) :
		_name(name),
		_type(type),
		_help(help),
		_func(func)
	{
	}

	void
	operator()() const
	{
		_func(_name);
	}

	const std::string&
	name() const
	{
		return _name;
	}

	Type
	type() const
	{
		return _type;
	}

	const std::string&
	help() const
	{
		return _help;
	}

private:
	std::string _name;
	Type _type{Type::invalid};
	std::string _help;
	std::function<void(const std::string&)> _func;
};

class JSONElemVector
{
public:
	JSONElemVector() = default;

	template <typename... Ts>
	JSONElemVector&
	add(Ts&&... args)
	{
		_elems.emplace_back(std::forward<Ts>(args)...);
		return *this;
	}

	void
	processParams();

	void
	print(std::ostream& os) const;

private:
	std::deque<JSONElem> _elems;
};
} // namespace detail
} // namespace options
} // namespace xolotl
