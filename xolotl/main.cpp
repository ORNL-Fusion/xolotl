#include <cstdlib>
#include <iostream>

#include <xolotl/interface/Interface.h>
#include <xolotl/options/CommandLineError.h>

int
main(int argc, const char* argv[])
{
	try {
		auto xi = xolotl::interface::makeXolotlInterface(argc, argv);
		if (xi) {
			xi->solveXolotl();
		}
	}
	catch (const xolotl::options::CommandLineError& e) {
		std::cout << e.what() << std::endl;
	}
	return EXIT_SUCCESS;
}
