#ifndef AMBIENT_P_ACTION_H
#define AMBIENT_P_ACTION_H
#include "ambient/p_profile.h"

namespace ambient
{
    class p_action {
    public:
        template <typename L, typename R>
        p_action(char op_code, const L& lhs, const R& rhs);
        const p_profile profile() const;
        std::pair<p_profile,p_profile> arguments;
        char op_code;
    };
}
#endif
