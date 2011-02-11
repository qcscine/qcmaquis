#ifndef AMBIENT_AUX_H
#define AMBIENT_AUX_H

#include "ambient/core/coherency.h"
#include <vector>

#define HASH_MAP_PARTITION_BIT_SIZE 8
#define HASH_MAP_PARTITION_SIZE 256 // 2^HAH_MAP_PARTITION_BIT_SIZE
#define HASH_MAP_VECTOR_RESERVATION 10

namespace ambient{

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

    template <class T>
    class hash_map {
    public:
        hash_map();
        void insert(unsigned int* hash, unsigned int hash_len, unsigned int id, T value, int shift = 0);
        T find(unsigned int* hash, unsigned int hash_len, unsigned int id, int shift = 0);
    private:
        std::vector< std::pair<hash_map*,std::vector<T>* > > content;
    };
   
}

#endif
