/*****************************************************************************
 *
 * MAQUIS Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_SIZEOF_H
#define MAQUIS_SIZEOF_H

namespace utils
{   
    template<class Iterator>
    std::size_t size_of(Iterator i1, Iterator i2)
    {
        std::size_t r = 0;
        for ( ; i1 != i2; ++i1)
            r += size_of(*i1);
        return r;
    }
}

#endif
