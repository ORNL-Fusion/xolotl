#pragma once

#include <boost/test/unit_test.hpp>

#include <Kokkos_Core.hpp>

class KokkosFixture
{
public:
	KokkosFixture() :
		_guard(boost::unit_test::framework::master_test_suite().argc,
			boost::unit_test::framework::master_test_suite().argv)
	{
	}

private:
	Kokkos::ScopeGuard _guard;
};
