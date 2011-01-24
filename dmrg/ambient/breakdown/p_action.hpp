#include "ambient/breakdown/p_action.h"
#include "ambient/breakdown/profiles.h"

namespace ambient { namespace breakdown {
    template <typename L, typename R>
    p_action::p_action(char op_code, const L& lhs, const R& rhs): arguments(get_profile<L>(lhs), get_profile<R>(rhs)), op_code(op_code) {
        scheduler::instance().push(this);
    };
    const p_profile p_action::profile() const { return p_profile(this, "action"); }
} }
