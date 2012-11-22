/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef NOOP_STORAGE_H
#define NOOP_STORAGE_H

class NoopStorage { };

class NoopStorageMaster { 
public:
    NoopStorage child() { return NoopStorage(); }
    void sync() { }
    void print_size() const { }
};

namespace storage
{
    template<class T> void store(T &, NoopStorage &) { }
    template<class T> void prefetch(T &, NoopStorage &) { }
    template<class T> void load(T &, NoopStorage &) { }
    void reset(NoopStorage &) { }
}

#endif
