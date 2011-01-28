#ifndef AMBIENT_AUX_H
#define AMBIENT_AUX_H

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

}

#endif
