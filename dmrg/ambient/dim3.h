#ifndef AMBIENT_DIM3_H
#define AMBIENT_DIM3_H

namespace ambient {

    class dim3
    {
    public:
        unsigned int x, y, z;
        dim3(unsigned int x = 1, unsigned int y = 1, unsigned int z = 1) : x(x), y(y), z(z) {}
        dim3& operator=(int value){
            x = y = z = value;
        }
        bool operator==(int value){
            return (x == value && y == value && z == value);
        }
    };

}
#endif
