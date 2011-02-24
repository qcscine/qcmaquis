#ifndef AMBIENT_AUX_H
#define AMBIENT_AUX_H

#include "ambient/core/layout.h"
#include <vector>

#define HASH_MAP_PARTITION_BIT_SIZE 8
#define HASH_MAP_PARTITION_SIZE 256 // 2^HAH_MAP_PARTITION_BIT_SIZE
#define HASH_MAP_VECTOR_RESERVATION 10
#define STACK_CONTENT_RESERVATION 10


namespace ambient{
namespace core{ class operation; }

    class dim3
    {
    public:
        unsigned int x, y, z;
        dim3(unsigned int x = 1, unsigned int y = 1, unsigned int z = 1) : x(x), y(y), z(z) {}
        dim3& operator=(int value){
            x = y = z = value;
        }
        dim3 operator*=(const dim3 & b){
            this->x *= b.x;
            this->y *= b.y;
            this->z *= b.z;
            return *this;
        }
        const dim3 operator*(const dim3 & b) const {
            return dim3(this->x, this->y, this->z) *= b;
        }
        bool operator==(int value){
            return (x == value && y == value && z == value);
        }
    };

    class hash_map {
    private: 
        hash_map();                             // constructor is private
        hash_map(hash_map const&);              // copy constructor is private
        hash_map& operator=(hash_map const&);   // assignment operator is private
    public:
        static hash_map& instance();
    public:
        unsigned int insert(unsigned int* hash, unsigned int hash_len, core::layout_table* value, int shift = 0);
        core::layout_table* find(unsigned int* hash, unsigned int hash_len, unsigned int id, int shift = 0) const;
    private:
        std::vector< std::pair<hash_map*,std::vector<core::layout_table*>* > > content;
    };

    class operation_stack {
    public:
        operation_stack();
        void push_back(std::pair<core::operation*,core::operation*> element);
        bool end_reached();
        std::pair<core::operation*,core::operation*>* pick();
    private:
        std::pair<core::operation*,core::operation*>* content;
        size_t write_iterator; 
        size_t read_iterator;
        size_t length;
        size_t reserved;
    };
}

#endif
