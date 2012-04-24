#ifndef AMBIENT_MODELS_OPERATION_SELECT
#define AMBIENT_MODELS_OPERATION_SELECT
#include <memory.h>
#include "ambient/ambient.h"

namespace ambient {
    template<int TT>
    int parseout(const char* sql, char** id);
    void ctxt_select(const char* sql);
    void ctxt_retain(const char* sql);
}

#endif
