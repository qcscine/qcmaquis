#ifndef MAQUIS_ALLOCATOR_HPP
#define MAQUIS_ALLOCATOR_HPP

#include <cstdlib>
#include <new>
#include <memory>

#if __cplusplus >= 201103L 
#define NOEXCEPT_SPEC noexcept
#else
#define NOEXCEPT_SPEC
#endif

namespace maquis {

template <typename T, unsigned int Size>
class static_allocator {
  public:
    using pointer = T *;
    using const_pointer = const T *;
    using reference = T &;
    using const_reference = const T &;
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template <typename U>
    struct rebind {
        using other = static_allocator<U, Size>;
    };

    static_allocator() NOEXCEPT_SPEC = default;

    static_allocator(static_allocator const& a) NOEXCEPT_SPEC { }

    template <typename U>
    static_allocator(static_allocator<U,Size> const& b) NOEXCEPT_SPEC { }

    pointer allocate(size_type n) {
        if (n > Size)
            throw std::bad_alloc();
        else
            return data;
    }

    void deallocate(pointer p, size_type n) NOEXCEPT_SPEC { }

    size_type max_size() const NOEXCEPT_SPEC {
        return Size;
    }

#if __cplusplus >= 201103L 
    template <typename C, class... Args>
    void construct(C* c, Args&&... args) {
        new ((void*)c) C(std::forward<Args>(args)...);
    }
#else
    void construct(pointer p, const_reference t) {
        new((void *)p) T(t);
    }
#endif

    template <typename C>
    void destroy(C* c) {
        c->~C();
    }

    bool operator == (static_allocator const& a2) const NOEXCEPT_SPEC {
        return true;
    }

    bool operator != (static_allocator const& a2) const NOEXCEPT_SPEC {
        return false;
    }

    template <typename U, unsigned int USize>
    bool operator == (static_allocator<U,USize> const& b) const NOEXCEPT_SPEC {
        return false;
    }

    template <typename U, unsigned int USize>
    bool operator != (static_allocator<U,USize> const& b) const NOEXCEPT_SPEC {
        return true;
    }

    private:
        T data[Size];
};

// Alignment must be a power of 2 !
template <typename T, unsigned int Alignment>
class aligned_allocator {
  public:
    using pointer = T *;
    using const_pointer = const T *;
    using reference = T &;
    using const_reference = const T &;
    using value_type = T;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    template <typename U>
    struct rebind {
        using other = aligned_allocator<U, Alignment>;
    };

    aligned_allocator() NOEXCEPT_SPEC = default;

    aligned_allocator(aligned_allocator const& a) NOEXCEPT_SPEC = default;

    template <typename U>
    aligned_allocator(aligned_allocator<U,Alignment> const& b) NOEXCEPT_SPEC {
    }

    pointer allocate(size_type n) {
        pointer p;
#ifdef _WIN32
        p = _aligned_malloc(n*sizeof(T), Alignment);
        if(p == 0)
            throw std::bad_alloc();
#else
        if(posix_memalign(reinterpret_cast<void**>(&p), Alignment, n * sizeof(T) ))
            throw std::bad_alloc();
#endif
        return p;
    }

    void deallocate(pointer p, size_type n) NOEXCEPT_SPEC {
        std::free(p);
    }

    size_type max_size() const NOEXCEPT_SPEC {
        std::allocator<T> a;
        return a.max_size();
    }

#if __cplusplus >= 201103L 
    template <typename C, class... Args>
    void construct(C* c, Args&&... args) {
        new ((void*)c) C(std::forward<Args>(args)...);
    }
#else
    void construct(pointer p, const_reference t) {
        new((void *)p) T(t);
    }
#endif

    template <typename C>
    void destroy(C* c) {
        c->~C();
    }

    bool operator == (aligned_allocator const& a2) const NOEXCEPT_SPEC {
        return true;
    }

    bool operator != (aligned_allocator const& a2) const NOEXCEPT_SPEC {
        return false;
    }

    template <typename U, unsigned int UAlignment>
    bool operator == (aligned_allocator<U,UAlignment> const& b) const NOEXCEPT_SPEC {
        return false;
    }

    template <typename U, unsigned int UAlignment>
    bool operator != (aligned_allocator<U,UAlignment> const& b) const NOEXCEPT_SPEC {
        return true;
    }
};

}

#undef NOEXPECT_SPEC

#endif //HPC12_ALIGNED_ALLOCATOR_HPP
