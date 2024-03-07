#include <xolotl/interface/Interface.h>
#include <cstdlib>

int
main(int argc, const char* argv[])
{
	xolotl::interface::makeXolotlInterface(argc, argv)->solveXolotl();
	return EXIT_SUCCESS;
}
