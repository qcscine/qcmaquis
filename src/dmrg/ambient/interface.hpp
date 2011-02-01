#ifndef AMBIENT_INTERFACE_H
#define AMBIENT_INTERFACE_H
#include <assert.h>
#include "ambient/core/operation.h"
#include "ambient/interface/p_profile.h"
#include "ambient/interface/charge.hpp"

namespace blas{ using namespace ambient;
#include "ambient/profiles.hpp"
#include "ambient/kernels.cpp"
} namespace ambient { using namespace blas;
#include "ambient/interface/core.hpp"
}

#endif
