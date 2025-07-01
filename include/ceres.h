#ifndef CERES_WRAPPER_H
#define CERES_WRAPPER_H

// Save any existing CERES_EXPORT definition
#pragma push_macro("CERES_EXPORT")
#pragma push_macro("CERES_NO_EXPORT")

// Undefine and redefine for static compilation
#ifdef CERES_EXPORT
#undef CERES_EXPORT
#endif
#ifdef CERES_NO_EXPORT  
#undef CERES_NO_EXPORT
#endif

// Define as empty for static compilation
#define CERES_EXPORT
#define CERES_NO_EXPORT

// WASM threading stubs - only provide what's absolutely needed
// Remove atomic stub since it conflicts with Eigen

// Include ceres headers
#include "ceres/version.h"
#include "ceres/ceres.h"

// Restore original definitions
#pragma pop_macro("CERES_NO_EXPORT")
#pragma pop_macro("CERES_EXPORT")

#endif // CERES_WRAPPER_H
