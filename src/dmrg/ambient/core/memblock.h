#ifndef AMBIENT_CORE_MEMBLOCK_H
#define AMBIENT_CORE_MEMBLOCK_H
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    class p_profile;

    class memblock {
    public:
        memblock(p_profile** p, int i, int j = 0);
        void set_memory(void* memory);
        void* element(int i, int j = 0);
        void* operator()(int i, int j = 0);
        void* item(int i, int j = 0);
        p_profile** profile;
        p_profile* get_profile();
        dim2 get_mem_dim();
        dim2 get_mem_t_dim();
        dim2 get_item_dim();
        void* header;
        void* data;
        size_t timestamp;
        bool available();
    private:
        template<typename T> operator T ();
    public:
        template<typename T> operator T* ()
        { return (T*)this->data;          }
        int owner;
        int i,j;
    };
}
#endif
