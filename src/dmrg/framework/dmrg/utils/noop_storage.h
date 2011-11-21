#ifndef NOOP_STORAGE_H
#define NOOP_STORAGE_H

class NoopStorage { };

class NoopStorageMaster { 
public:
    NoopStorage child() { return NoopStorage(); }

};

namespace storage
{
    template<class T> void store(T &, NoopStorage &) { }
    template<class T> void prefetch(T &, NoopStorage &) { }
    template<class T> void load(T &, NoopStorage &) { }
}

#endif

