/**
 * @file test_api.c
 * @brief Bootstraps the test framework for the CS248 library.
 * @author - SCEC
 * @version 1.0
 *
 * Tests the cS248 library by loading it and executing the code as
 * UCVM would do it.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "cs248.h"

/**
 * Initializes and runs the test program. Tests link against the
 * static version of the library to prevent any dynamic loading
 * issues.
 *
 * @param argc The number of arguments.
 * @param argv The argument strings.
 * @return A zero value indicating success.
 */
int main(int argc, const char* argv[]) {

	// Declare the structures.
	cs248_point_t pt;
	cs248_properties_t ret;

        // Initialize the model.
        char *envstr=getenv("UCVM_INSTALL_PATH");
        if(envstr != NULL) {
           assert(cs248_init(envstr, "cs248") == 0);
           } else {
             assert(cs248_init("..", "cs248") == 0);
        }

	printf("Loaded the model successfully.\n");

	// Query a point.
	pt.longitude = -120;
	pt.latitude = 37.2;
	pt.depth = 0;

	cs248_query(&pt, &ret, 1);

	assert(ret.vs > 0);
	assert(ret.vp > 0);
	assert(ret.rho > 0);

	printf("Query was successful.\n");

	// Close the model.
	assert(cs248_finalize() == 0);

	printf("Model closed successfully.\n");

	printf("\nALL CS248 TESTS PASSED\n");

	return 0;
}
