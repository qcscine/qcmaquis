/*
 * Ambient, License - Version 1.0 - May 3rd, 2012
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef AMBIENT_MEMORY_ALLOCATOR
#define AMBIENT_MEMORY_ALLOCATOR

namespace ambient {

    template <class T>
    class default_allocator {
    public:
        static void* alloc(pool::descriptor& spec){
            return ambient::pool::malloc(spec);
        }
        static void* calloc(pool::descriptor& spec){
            void* m = ambient::pool::malloc(spec);
            memset(m, 0, spec.extent);
            return m;
        }
        static void free(void* ptr, pool::descriptor& spec){
            ambient::pool::free(ptr, spec);
        }
        // std:: compatibility aliases //
        static T* allocate(std::size_t n){ 
        }
        static void deallocate(T* ptr, std::size_t n){ 
        }
        template <class U> struct rebind {
            typedef default_allocator<U> other; 
        };
        typedef T value_type;
    };

    template <class T>
    class constrained_allocator {
    public:
        static void* alloc(pool::descriptor& spec){
            if(spec.mmap == NULL) return ambient::pool::malloc(spec);
            return ambient::pool::malloc<outofcore>(spec);
        }
        static void* calloc(pool::descriptor& spec){
            void* m = alloc(spec);
            memset(m, 0, spec.extent);
            return m;
        }
        static void free(void* ptr, pool::descriptor& spec){
        } // lets leek for now // do we ?

        // std:: compatibility aliases //
        static T* allocate(std::size_t n){ 
        }
        static void deallocate(T* ptr, std::size_t n){ 
        }
        template <class U> struct rebind {
            typedef constrained_allocator<U> other; 
        };
        typedef T value_type;
    };

    template <class T>
    class bulk_allocator {
    public:
        static T* allocate(std::size_t n){
            return (T*)ambient::pool::malloc<bulk>(n*sizeof(T));
        }
        static void deallocate(T* p, std::size_t n){}
        template <class U> struct rebind { 
            typedef bulk_allocator<U> other; 
        };
        typedef T value_type;
    };

}

#endif
