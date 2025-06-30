// Test including ceres headers
#include "ceres/version.h"

// Simple function that uses ceres header
extern "C" int test_ceres() {
    return CERES_VERSION_MAJOR;
}
