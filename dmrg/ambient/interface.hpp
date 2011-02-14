#ifndef AMBIENT_INTERFACE_H
#define AMBIENT_INTERFACE_H
#include <assert.h>
#include "ambient/core/operation.h"
#include "ambient/core/p_profile.h"
#include "ambient/core/select.h"

namespace blas{ using namespace ambient;
typedef ambient::p_profile   void_pt;
typedef ambient::p_profile_s void_spt;
#include "ambient/kernels.cpp"
#include "ambient/interface/profiles.hpp"
} namespace ambient { using namespace blas;
#include "ambient/interface/core.hpp"
}

#endif
