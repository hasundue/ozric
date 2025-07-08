#pragma once

#ifdef CERES_EXPORT
#undef CERES_EXPORT
#endif

#ifdef CERES_NO_EXPORT
#undef CERES_NO_EXPORT
#endif

#define CERES_EXPORT
#define CERES_NO_EXPORT

#include <ceres/ceres.h>
