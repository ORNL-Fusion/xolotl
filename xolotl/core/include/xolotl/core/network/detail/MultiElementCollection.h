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
template <typename T>
struct TypeTag
{
	using Type = T;
};

template <typename TElem>
class ElementSet
{
public:
	using ElementType = TElem;
	using IndexType = ::xolotl::IdType;

	ElementSet() = default;

	ElementSet(Kokkos::View<ElementType*> elemView) :
		_elems(elemView),
		_numElems(elemView.size())
	{
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return _elems.required_allocation_size(_numElems);
	}

	IndexType
	getNumberOfElements() const noexcept
	{
		return _numElems;
	}

	Kokkos::View<ElementType*>
	getView() const
	{
		return _elems;
	}

	void
	setView(Kokkos::View<ElementType*> view)
	{
		_elems = view;
		_numElems = view.size();
	}

	template <typename F>
	KOKKOS_INLINE_FUNCTION
	void
	apply(const F& func, const IndexType i) const
	{
		func(_elems(i));
	}

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

template <std::size_t NumElementTypes, typename... TElems>
class ElementSetMixinChain
{
	template <std::size_t, typename...>
	friend class ElementSetMixinChain;

public:
	using IndexType = ::xolotl::IdType;

	ElementSetMixinChain() = default;

	ElementSetMixinChain(IndexType indexBegin) : _indexBegin(indexBegin)
	{
	}

	void updateIndices(IndexType)
	{
	}

	void
	getView() const
	{
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return 0;
	}

	IndexType
	getNumberOfElements() const noexcept
	{
		return 0;
	}

	void
	getElementBeginIndices(
		Kokkos::Array<IndexType, NumElementTypes + 1>& ids) const
	{
		ids[NumElementTypes] = _indexBegin;
	}

	template <typename F>
	void
	invoke(const F&) const noexcept
	{
	}

	template <typename F>
	KOKKOS_INLINE_FUNCTION
	void
	apply(const F&, const IndexType) const noexcept
	{
	}

	template <typename F, typename T>
	KOKKOS_INLINE_FUNCTION
	void
	reduce(const F&, const IndexType, T&) const
	{
	}

private:
	KOKKOS_INLINE_FUNCTION
	IndexType
	getIndexBegin() const noexcept
	{
		return _indexBegin;
	}

private:
	IndexType _indexBegin{};
};

template <std::size_t NumElementTypes, typename TElem, typename... TOtherElems>
class ElementSetMixinChain<NumElementTypes, TElem, TOtherElems...> :
	public ElementSet<TElem>,
	ElementSetMixinChain<NumElementTypes, TOtherElems...>
{
	template <std::size_t, typename...>
	friend class ElementSetMixinChain;

	static constexpr std::size_t
	getElementTypeIndex() noexcept
	{
		return NumElementTypes - (sizeof...(TOtherElems) + 1);
	}

public:
	using IndexType = ::xolotl::IdType;
	using Head = ElementSet<TElem>;
	using Tail = ElementSetMixinChain<NumElementTypes, TOtherElems...>;

	ElementSetMixinChain() = default;

	template <typename THeadView, typename... TTailViews>
	ElementSetMixinChain(
		IndexType indexBegin, THeadView view, TTailViews... views) :
		Head(view),
		Tail(indexBegin + static_cast<IndexType>(view.size()), views...),
		_indexBegin(indexBegin)
	{
	}

	template <typename... TViews>
	ElementSetMixinChain(TViews... views) :
		ElementSetMixinChain(static_cast<IndexType>(0), views...)
	{
	}

	void
	updateIndices(IndexType indexBegin)
	{
		_indexBegin = indexBegin;
		Tail::updateIndices(indexBegin + Head::getNumberOfElements());
	}

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

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return Head::getDeviceMemorySize() + Tail::getDeviceMemorySize();
	}

	IndexType
	getNumberOfElements() const noexcept
	{
		return Head::getNumberOfElements() + Tail::getNumberOfElements();
	}

	void
	getElementBeginIndices(
		Kokkos::Array<IndexType, NumElementTypes + 1>& ids) const
	{
		ids[getElementTypeIndex()] = _indexBegin;
		Tail::getElementBeginIndices(ids);
	}

	template <typename F>
	void
	invoke(const F& func) const
	{
		func(getElementTypeIndex(), Head::getNumberOfElements(),
			TypeTag<TElem>{});
		Tail::invoke(func);
	}

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
	KOKKOS_INLINE_FUNCTION
	IndexType
	getIndexBegin() const noexcept
	{
		return _indexBegin;
	}

private:
	IndexType _indexBegin{};
};

template <typename TElementTypeList>
struct ElementSetChainHelper;

template <typename... TElems>
struct ElementSetChainHelper<std::tuple<TElems...>>
{
	using Type = ElementSetMixinChain<sizeof...(TElems), TElems...>;
};

template <typename TElementTypeList>
using ElementSetChain = typename ElementSetChainHelper<TElementTypeList>::Type;

template <typename... TElems>
class MultiElementCollection
{
	static constexpr std::size_t numElementTypes = sizeof...(TElems);

public:
	using IndexType = ::xolotl::IdType;

	MultiElementCollection() = default;

	template <typename... TViews>
	MultiElementCollection(TViews... views) :
		_chain(views...),
		_numElems(_chain.getNumberOfElements())
	{
		static_assert(sizeof...(TViews) == numElementTypes,
			"Construction from views requires the number of views to match the "
			"number of element types in the MultiElementCollection");
	}

	std::uint64_t
	getDeviceMemorySize() const noexcept
	{
		return _chain.getDeviceMemorySize();
	}

	IndexType
	getNumberOfElements() const noexcept
	{
		return _numElems;
	}

	template <typename TElem>
	Kokkos::View<TElem*>
	getView() const
	{
		return _chain.template getView<TElem>();
	}

	template <typename TElem>
	void
	setView(Kokkos::View<TElem*> view)
	{
		_chain.setView(view);
		_numElems = _chain.getNumberOfElements();
	}

	template <typename F>
	void
	forEach(const F& func)
	{
		auto chain = _chain;
		Kokkos::parallel_for(
			_numElems,
			DEVICE_LAMBDA(const IndexType i) { chain.apply(func, i); });
	}

	template <typename TElem, typename F>
	void
	forEachOn(const F& func)
	{
		auto view = getView<TElem>();
		Kokkos::parallel_for(
			view.size(), DEVICE_LAMBDA(const IndexType i) { func(view[i]); });
	}

	template <typename F>
	void
	forEachType(const F& func)
	{
		_chain.invoke(func);
	}

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

	decltype(auto)
	getChain() const
	{
		return _chain;
	}

	Kokkos::Array<IndexType, numElementTypes + 1>
	getElementBeginIndices() const
	{
		Kokkos::Array<IndexType, numElementTypes + 1> ret;
		_chain.getElementBeginIndices(ret);
		return ret;
	}

private:
	ElementSetChain<std::tuple<TElems...>> _chain;
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
