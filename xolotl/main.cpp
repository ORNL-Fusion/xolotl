#include <xolotl/interface/Interface.h>

int
main(int argc, char* argv[])
{
	xolotl::interface::XolotlInterface{argc, argv}.solveXolotl();
	return EXIT_SUCCESS;
}
