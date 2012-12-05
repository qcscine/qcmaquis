namespace ambient { namespace models { namespace velvet {

    inline memspec::memspec(dim2 dim, size_t ts) : dim(dim), size(dim.square()*ts) { }

} } }
