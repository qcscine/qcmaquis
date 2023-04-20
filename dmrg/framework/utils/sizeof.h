/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
