#ifndef AMBIENT_CORE_SELECT_H
#define AMBIENT_CORE_SELECT_H

#include <memory.h>
#include "ambient/ambient.h"

namespace ambient {
    int parseout_id(const char* sql, char** id);
    void scope_select(const char* sql);
    void scope_retain(const char* sql);
}

#endif
