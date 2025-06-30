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

// Include ceres headers (use full path to avoid recursion)
#include "/home/hasundue/.cache/zig/p/N-V-__8AAAwZyQAma1OByLx2AKj3upoY7FxoS7bQ1LBLR2m3/include/ceres/version.h"
#include "/home/hasundue/.cache/zig/p/N-V-__8AAAwZyQAma1OByLx2AKj3upoY7FxoS7bQ1LBLR2m3/include/ceres/ceres.h"

// Restore original definitions
#pragma pop_macro("CERES_NO_EXPORT")
#pragma pop_macro("CERES_EXPORT")

#endif // CERES_WRAPPER_H