#ifndef AMBIENT_CORE_WORKGROUP_H
#define AMBIENT_CORE_WORKGROUP_H
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    class p_profile;

    class workgroup {
    public:
        workgroup(p_profile** p, int i, int j = 0, int k = 0);
        void set_memory(void* memory);
        void* element(int i, int j = 0, int k = 0);
        void* operator()(int i, int j = 0, int k = 0);
        void* item(int i, int j = 0, int k = 0);
        p_profile** profile;
        p_profile* get_profile();
        dim3 get_group_dim();
        dim3 get_group_t_dim();
        dim3 get_item_dim();
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
        int i, j, k;
    };
}
#endif
