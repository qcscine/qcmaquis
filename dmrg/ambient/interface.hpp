#ifndef AMBIENT_INTERFACE_H
#define AMBIENT_INTERFACE_H
#include <assert.h>
#include "ambient/core/operation.h"
#include "ambient/core/p_profile.h"
#include "ambient/core/select.h"

namespace blas{ using namespace ambient;
#include "ambient/interface/profiles.hpp"
#include "ambient/kernels.cpp"
} namespace ambient { using namespace blas;
#include "ambient/interface/core.hpp"
}

#endif
