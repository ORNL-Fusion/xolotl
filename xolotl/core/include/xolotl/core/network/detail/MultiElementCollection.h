#pragma once

#include <tuple>
#include <type_traits>

#include <Kokkos_Core.hpp>

#if defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#define DEVICE_LAMBDA [=] __device__
#else
#define DEVICE_LAMBDA [=]
#endif

#include <xolotl/core/network/detail/TupleUtility.h>

namespace xolotl
{
namespace core
{
namespace network
{
namespace detail
{
/**
 * @brief Wrap a Kokkos::View of a single element type and provide interface
 * for ElementSetMixinChain
 */
template <typename TElem>
class ElementSetMixin
{
public:
	using ElementType = TElem;
	using IndexType = ::xolotl::IdType;

	/**
	 * @brief Tag type for ElementType used in the invoke() function
	 */
	struct ElemTypeTag
	{
		using Type = ElementType;
	};

	/**
	 * @brief Default to empty set
	 */
	ElementSetMixin() = default;

	/**
	 * @brief Construct with populated view
	 */
	ElementSetMixin(Kokkos::View<ElementType*> elemView) :
		_elems(elemView),
		_numElems(elemView.size())
	{
	}

	/**
	 * @brief Get required device allocation size
	 */
	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return _elems.required_allocation_size(_numElems);
	}

	/**
	 * @brief Get size of current set
	 */
	IndexType
	getNumberOfElements() const noexcept
	{
		return _numElems;
	}

	/**
	 * @brief Get wrapped view
	 */
	Kokkos::View<ElementType*>
	getView() const
	{
		return _elems;
	}

	/**
	 * @brief Replace wrapped view
	 */
	void
	setView(Kokkos::View<ElementType*> view)
	{
		_elems = view;
		_numElems = view.size();
	}

	/**
	 * @brief Invoke the given function with the element type tag and index
	 *
	 * This is for invoking a function for each element type rather than for
	 * each element. This is part of the implementation of
	 * MultiElementCollection::forEachType().
	 */
	template <typename F>
	void
	invoke(const F& func, std::size_t elementTypeIndex) const
	{
		func(elementTypeIndex, getNumberOfElements(), ElemTypeTag{});
	}

	/**
	 * @brief Invoke the given function with the i-th element of this set
	 */
	template <typename F>
	KOKKOS_INLINE_FUNCTION
	void
	apply(const F& func, const IndexType i) const
	{
		func(_elems(i));
	}

	/**
	 * @brief Invoke the given (reduction) function with the i-th element of
	 * this set
	 */
	template <typename F, typename T>
	KOKKOS_INLINE_FUNCTION
	void
	reduce(const F& func, const IndexType i, T& local) const
	{
		func(_elems(i), local);
	}

private:
	Kokkos::View<ElementType*> _elems;
	IndexType _numElems{};
};

/**
 * @brief The un-specialized version is used as end-point for recursively
 * constructed inheritance chain
 *
 * This class provides the same interface as its specialized counterpart, but
 * these functions do nothing since this class does not own any data.
 *
 * Besides stubbing the interface, this class has a begin index which is
 * expected to be one-past-the-last index for the entire MultiElementCollection,
 * since this class is used only at the end of the chain
 *
 * @tparam NumElementTypes Total number of element types used in the chain
 * @tparam TElems Parameter pack of element types (expected to be empty for
 * default instantiation)
 */
template <std::size_t NumElementTypes, typename... TElems>
class ElementSetMixinChain
{
	template <std::size_t, typename...>
	friend class ElementSetMixinChain;

public:
	using IndexType = ::xolotl::IdType;

	/**
	 * @brief Default sets begin index to 0
	 */
	ElementSetMixinChain() = default;

	/**
	 * @brief Construct with begin index
	 */
	ElementSetMixinChain(IndexType indexBegin) : _indexBegin(indexBegin)
	{
	}

	/**
	 * @brief Does nothing
	 */
	void
	updateIndices(IndexType)
	{
	}

	/**
	 * @brief Does nothing
	 */
	void
	getView() const
	{
	}

	/**
	 * @brief Returns zero
	 */
	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return 0;
	}

	/**
	 * @brief Returns zero
	 */
	IndexType
	getNumberOfElements() const noexcept
	{
		return 0;
	}

	/**
	 * @brief The begin index for this version is expected to be
	 * one-past-the-last for the entire collection since an instance of this
	 * version is at the end of the chain
	 */
	void
	getElementBeginIndices(
		Kokkos::Array<IndexType, NumElementTypes + 1>& ids) const
	{
		ids[NumElementTypes] = _indexBegin;
	}

	/**
	 * @brief Does nothing
	 */
	template <typename F>
	void
	invoke(const F&) const noexcept
	{
	}

	/**
	 * @brief Does nothing
	 */
	template <typename F>
	KOKKOS_INLINE_FUNCTION
	void
	apply(const F&, const IndexType) const noexcept
	{
	}

	/**
	 * @brief Does nothing
	 */
	template <typename F, typename T>
	KOKKOS_INLINE_FUNCTION
	void
	reduce(const F&, const IndexType, T&) const
	{
	}

private:
	/**
	 * @brief Get begin index
	 */
	KOKKOS_INLINE_FUNCTION
	IndexType
	getIndexBegin() const noexcept
	{
		return _indexBegin;
	}

private:
	IndexType _indexBegin{};
};

/**
 * @brief The specialized version is instantiated recursively for each element
 * type and chained together by inheritance.
 *
 * This class handles the bulk of the implementation machinery using recursion
 * logic. Given a list of element types, an instantiation of this class takes
 * the first element type to instantiate an ElementSetMixin to own the data for
 * that type. The rest of the types in the list are passed up to another
 * instantiation of this class that we inherit from. This constructs a custom
 * inheritance hierarchy at compile time. When the "tail" is empty, the chain
 * will be capped by the default instantation which does not manage any
 * elements.
 *
 * Here's an example of an inheritance graph given the three types `P, D, R`.
 *
 * ```
 *                                 ElementSetMixin<R>    ElementSetMixinChain<3>
 *                                       \                  /
 *                                        \    ____________/
 *                                         \  /
 *                                          \/
 *                 ElementSetMixin<D>    ElementSetMixinChain<3, R>
 *                       \                  /
 *                        \    ____________/
 *                         \  /
 *                          \/
 * ElementSetMixin<P>    ElementSetMixinChain<3, D, R>
 *       \                  /
 *        \    ____________/
 *         \  /
 *          \/
 *        ElementSetMixinChain<3, P, D, R>
 * ```
 *
 * @tparam NumElementTypes Total number of element types in the
 * MultiElementCollection
 * @tparam TElem The element type used in the current instantiation's
 * ElementSetMixin
 * @tparam TOtherElems Zero or more other element types to be used to
 * instantiate other links above in the chain
 */
template <std::size_t NumElementTypes, typename TElem, typename... TOtherElems>
class ElementSetMixinChain<NumElementTypes, TElem, TOtherElems...> :
	public ElementSetMixin<TElem>,
	ElementSetMixinChain<NumElementTypes, TOtherElems...>
{
	template <std::size_t, typename...>
	friend class ElementSetMixinChain;

	/**
	 * @brief Get the index (position) of current element type within the list
	 * of element types provided to MultiElementCollection
	 */
	static constexpr std::size_t
	getElementTypeIndex() noexcept
	{
		return NumElementTypes - (sizeof...(TOtherElems) + 1);
	}

public:
	using IndexType = ::xolotl::IdType;
	/**
	 * @brief Short-hand referring to the set managing the first element
	 * type in the list
	 */
	using Head = ElementSetMixin<TElem>;
	/**
	 * @brief Short-hand referring to the rest of the chain above
	 */
	using Tail = ElementSetMixinChain<NumElementTypes, TOtherElems...>;

	/**
	 * @brief Default to empty
	 */
	ElementSetMixinChain() = default;

	/**
	 * @brief Construct with begin index, view to wrap, and other views to pass
	 * up the chain
	 */
	template <typename THeadView, typename... TTailViews>
	ElementSetMixinChain(
		IndexType indexBegin, THeadView view, TTailViews... views) :
		Head(view),
		Tail(indexBegin + static_cast<IndexType>(view.size()), views...),
		_indexBegin(indexBegin)
	{
	}

	/**
	 * @brief Construct from list of views to wrap
	 *
	 * This should be used only for the first instantiation of this class in the
	 * chain. This initiates the construction of the chain
	 */
	template <typename... TViews>
	ElementSetMixinChain(TViews... views) :
		ElementSetMixinChain(static_cast<IndexType>(0), views...)
	{
	}

	/**
	 * @brief Recursively set begin indices for each set in the chain based on
	 * the number of elements of each type
	 */
	void
	updateIndices(IndexType indexBegin)
	{
		_indexBegin = indexBegin;
		Tail::updateIndices(indexBegin + Head::getNumberOfElements());
	}

	/**
	 * @brief get the view owning the elements of type TGetElem
	 */
	template <typename TGetElem>
	decltype(auto)
	getView() const
	{
		if constexpr (std::is_same_v<TGetElem, TElem>) {
			return Head::getView();
		}
		else {
			return Tail::template getView<TGetElem>();
		}
	}

	template <typename TSetElem>
	void
	setView(Kokkos::View<TSetElem*> view)
	{
		if constexpr (std::is_same_v<TSetElem, TElem>) {
			Head::setView(view);
			updateIndices(_indexBegin);
		}
		else {
			Tail::setView(view);
		}
	}

	/**
	 * @brief Get required device allocation size
	 */
	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return Head::getDeviceMemorySize() + Tail::getDeviceMemorySize();
	}

	/**
	 * @brief Get the total number of elements managed by the current
	 * instantiation and those above in the chain
	 */
	IndexType
	getNumberOfElements() const noexcept
	{
		return Head::getNumberOfElements() + Tail::getNumberOfElements();
	}

	/**
	 * @brief Provide begin index for current instantiation's elements and pass
	 * on for the rest of the chain to do the same
	 */
	void
	getElementBeginIndices(
		Kokkos::Array<IndexType, NumElementTypes + 1>& ids) const
	{
		ids[getElementTypeIndex()] = _indexBegin;
		Tail::getElementBeginIndices(ids);
	}

	/**
	 * @brief Invoke the given function for the current element type, then pass
	 * up the chain to be invoked for the remaining element types
	 */
	template <typename F>
	void
	invoke(const F& func) const
	{
		Head::invoke(func, getElementTypeIndex());
		Tail::invoke(func);
	}

	/**
	 * @brief Invoke the given function with the appropriate element from the
	 * collection
	 *
	 * If the index belongs to the current instantiation's set of elements, we
	 * will apply the function to the corresponding element. Otherwise, we will
	 * pass up the chain.
	 */
	template <typename F>
	KOKKOS_INLINE_FUNCTION
	void
	apply(const F& func, const IndexType i) const
	{
		if (i < Tail::getIndexBegin()) {
			Head::apply(func, i - _indexBegin);
		}
		else {
			Tail::apply(func, i);
		}
	}

	/**
	 * @brief Invoke the given reduction function with the appropriate element
	 * from the collection
	 *
	 * If the index belongs to the current instantiation's set of elements, we
	 * will apply the function to the corresponding element. Otherwise, we will
	 * pass up the chain.
	 */
	template <typename F, typename T>
	KOKKOS_INLINE_FUNCTION
	void
	reduce(const F& func, const IndexType i, T& local) const
	{
		if (i < Tail::getIndexBegin()) {
			Head::reduce(func, i - _indexBegin, local);
		}
		else {
			Tail::reduce(func, i, local);
		}
	}

private:
	/**
	 * @brief Get the begin index for the current element type
	 */
	KOKKOS_INLINE_FUNCTION
	IndexType
	getIndexBegin() const noexcept
	{
		return _indexBegin;
	}

private:
	IndexType _indexBegin{};
};

/**
 * @brief Linear collection of multiple element types that behaves as single
 * container
 *
 * Intuitively, views of elements of different types are connected end-to-end so
 * that a single index can be used (between zero and the total number of
 * elements) to access any element in the collection regardless of type.
 *
 * The motivation for this kind of container is to perform a similar operation
 * on multiple element sets with only one GPU kernel.
 *
 * The bulk of the logic is implemented through
 * ElementSetMixinChain< NumElementTypes, TElem, TOtherElems... >.
 */
template <typename... TElems>
class MultiElementCollection
{
	static constexpr std::size_t numElementTypes = sizeof...(TElems);

public:
	using IndexType = ::xolotl::IdType;

	/**
	 * @brief Default construct to empty collection
	 */
	MultiElementCollection() = default;

	/**
	 * @brief Construct with views of elements to be wrapped
	 *
	 * @note The views must of the same element types managed by the collection
	 * and in the same order as those element type template arguments
	 */
	template <typename... TViews>
	MultiElementCollection(TViews... views) :
		_chain(views...),
		_numElems(_chain.getNumberOfElements())
	{
		static_assert(sizeof...(TViews) == numElementTypes,
			"Construction from views requires the number of views to match the "
			"number of element types in the MultiElementCollection");
	}

	/**
	 * @brief Get required device allocation size
	 */
	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return _chain.getDeviceMemorySize();
	}

	/**
	 * @brief Get the total number of elements in the collection
	 */
	IndexType
	getNumberOfElements() const noexcept
	{
		return _numElems;
	}

	/**
	 * @brief Get the view owning the elements of type TElem
	 */
	template <typename TElem>
	Kokkos::View<TElem*>
	getView() const
	{
		return _chain.template getView<TElem>();
	}

	/**
	 * @brief Provide a view of elements of type TElem
	 */
	template <typename TElem>
	void
	setView(Kokkos::View<TElem*> view)
	{
		_chain.setView(view);
		_numElems = _chain.getNumberOfElements();
	}

	/**
	 * @brief Perform a Kokkos parallel_for on all the elements in the
	 * collection
	 *
	 * The callable should be of the form `void f(ElemType&& elem)` and
	 * templated on the type of the element parameter.
	 * This can be a generic lambda or a functor with a template call operator.
	 */
	template <typename F>
	void
	forEach(const F& func)
	{
		auto chain = _chain;
		Kokkos::parallel_for(
			_numElems,
			DEVICE_LAMBDA(const IndexType i) { chain.apply(func, i); });
	}

	template <typename F>
	void
	forEach(const std::string& label, const F& func)
	{
		auto chain = _chain;
		Kokkos::parallel_for(
			label, _numElems,
			DEVICE_LAMBDA(const IndexType i) { chain.apply(func, i); });
	}

	/**
	 * @brief Perform a Kokkos parallel_for on all the elements of a single type
	 */
	template <typename TElem, typename F>
	void
	forEachOn(const F& func)
	{
		auto view = getView<TElem>();
		Kokkos::parallel_for(
			view.size(), DEVICE_LAMBDA(const IndexType i) { func(view[i]); });
	}

	template <typename TElem, typename F>
	void
	forEachOn(const std::string& label, const F& func)
	{
		auto view = getView<TElem>();
		Kokkos::parallel_for(
			label, view.size(),
			DEVICE_LAMBDA(const IndexType i) { func(view[i]); });
	}

	/**
	 * @brief Apply a function for each distinct element type
	 */
	template <typename F>
	void
	forEachType(const F& func)
	{
		_chain.invoke(func);
	}

	/**
	 * @brief Perform a Kokkos parallel_reduce on all the elements in the
	 * collection
	 *
	 * The callable should be of the form `void f(ElemType&& elem, T& local)`
	 * and templated on the type of the element parameter.
	 * This can be a generic lambda or a functor with a template call operator.
	 *
	 * @tparam F Type of callable
	 * @tparam T Type of reduction variable
	 */
	template <typename F, typename T>
	void
	reduce(const F& func, T& out)
	{
		auto chain = _chain;
		Kokkos::parallel_reduce(
			_numElems,
			DEVICE_LAMBDA(
				const IndexType i, T& local) { chain.reduce(func, i, local); },
			out);
	}

	template <typename F, typename T>
	void
	reduce(const std::string& label, const F& func, T& out)
	{
		auto chain = _chain;
		Kokkos::parallel_reduce(
			label, _numElems,
			DEVICE_LAMBDA(
				const IndexType i, T& local) { chain.reduce(func, i, local); },
			out);
	}

	/**
	 * @brief Perform a Kokkos parallel_reduce on all the elements of a single
	 * type
	 */
	template <typename TElem, typename F, typename T>
	void
	reduceOn(const F& func, T& out)
	{
		auto view = getView<TElem>();
		Kokkos::parallel_reduce(
			view.size(),
			DEVICE_LAMBDA(
				const IndexType i, T& local) { func(view[i], local); },
			out);
	}

	template <typename TElem, typename F, typename T>
	void
	reduceOn(const std::string& label, const F& func, T& out)
	{
		auto view = getView<TElem>();
		Kokkos::parallel_reduce(
			label, view.size(),
			DEVICE_LAMBDA(
				const IndexType i, T& local) { func(view[i], local); },
			out);
	}

	/**
	 * @brief Get the ElementSetMixinChain that implements this collection
	 */
	decltype(auto)
	getChain() const
	{
		return _chain;
	}

	/**
	 * @brief Get the set of begin indices for each element view as if the views
	 * were appended one to another
	 */
	Kokkos::Array<IndexType, numElementTypes + 1>
	getElementBeginIndices() const
	{
		Kokkos::Array<IndexType, numElementTypes + 1> ret;
		_chain.getElementBeginIndices(ret);
		return ret;
	}

private:
	ElementSetMixinChain<sizeof...(TElems), TElems...> _chain;
	std::size_t _numElems{};
};

template <typename... TElems>
class MultiElementCollection<std::tuple<TElems...>> :
	public MultiElementCollection<TElems...>
{
	using Superclass = MultiElementCollection<TElems...>;

	using Superclass::Superclass;
};
} // namespace detail
} // namespace network
} // namespace core
} // namespace xolotl
