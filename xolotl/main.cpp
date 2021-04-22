#include <xolotl/interface/Interface.h>

int
main(int argc, const char* argv[])
{
	xolotl::interface::XolotlInterface{argc, argv}.solveXolotl();
	return EXIT_SUCCESS;
}
