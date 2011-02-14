#ifndef AMBIENT_INTERFACE_SELECT_H
#define AMBIENT_INTERFACE_SELECT_H

#include <memory.h>
#include "ambient/ambient.h"

namespace ambient {

    int parseout_id(const char* sql, char** id);
    void select(const char* sql);

}

#endif
