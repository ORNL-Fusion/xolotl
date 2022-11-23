#pragma once

#include <boost/test/unit_test.hpp>

#include <Kokkos_Core.hpp>

class KokkosFixture
{
public:
	KokkosFixture() :
		_guard((!Kokkos::is_initialized() && !Kokkos::is_finalized()) ?
				std::make_unique<Kokkos::ScopeGuard>(
					boost::unit_test::framework::master_test_suite().argc,
					boost::unit_test::framework::master_test_suite().argv) :
				nullptr)
	{
	}

private:
	std::unique_ptr<Kokkos::ScopeGuard> _guard;
};
