#ifndef TESTS_UTILS_MPIFIXTURE_H
#define TESTS_UTILS_MPIFIXTURE_H

// Do not put this in a namespace - it causes macros like 
// BOOST_GLOBAL_FIXTURE to break.

struct MPIFixture {
    MPIFixture(void) {
        BOOST_TEST_MESSAGE("Initializing MPI");
        auto& mts = boost::unit_test::framework::master_test_suite();
        MPI_Init(&mts.argc, &mts.argv);
        BOOST_TEST_MESSAGE("Done initializing MPI");
    }

    ~MPIFixture(void) {
        BOOST_TEST_MESSAGE("Finalizing MPI");
        MPI_Finalize();
        BOOST_TEST_MESSAGE("Done finalizing MPI");
    }
};

#endif // TESTS_UTILS_MPIFIXTURE_H
