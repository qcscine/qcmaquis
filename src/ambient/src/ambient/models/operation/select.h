#ifndef AMBIENT_CORE_SELECT_H
#define AMBIENT_CORE_SELECT_H

#include <memory.h>
#include "ambient/ambient.h"

namespace ambient {
    template<int TT>
    int parseout(const char* sql, char** id);
    void ctxt_select(const char* sql);
    void ctxt_retain(const char* sql);
}

#endif
