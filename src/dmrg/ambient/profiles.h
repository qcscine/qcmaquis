#ifndef AMBIENT_PROFILES_H
#define AMBIENT_PROFILES_H
#include "ambient/p_profile.h"

namespace ambient {

    template<typename T>
    const p_profile get_profile(const T& obj){ return obj.profile(); }
    template<>
    const p_profile get_profile<>(const int& obj){ return p_profile(&obj, "int"); }
    template<>
    const p_profile get_profile<>(const double& obj){ return p_profile(&obj, "double"); }

}
#endif
