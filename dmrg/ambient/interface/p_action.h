#ifndef AMBIENT_P_ACTION_H
#define AMBIENT_P_ACTION_H

namespace ambient {

    class p_profile;
    class p_action {
    public:
        template <typename L, typename R>
        p_action(char op_code, const L* lhs, const R* rhs);
        std::pair<const p_profile*,const p_profile*> arguments;
        char op_code;
    };
}
#endif
