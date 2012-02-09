#ifndef AMBIENT_DIM2_H
#define AMBIENT_DIM2_H
namespace ambient{

    class dim2 {
    public:
        size_t x, y;
        dim2(size_t x = 1, size_t y = 1) 
        : x(x), y(y) 
        {
        }

        dim2& operator=(int value){
            this->x = this->y = value;
            return *this;
        }

        dim2& operator*=(const dim2 & b){
            this->x *= b.x;
            this->y *= b.y;
            return *this;
        }

        size_t operator*(const dim2 & b) const { // multiplication of all components
            return this->x * b.x *
                   this->y * b.y ;
        }

        bool operator==(int value) const {
            return (x == value && y == value);
        }

        bool operator==(dim2 m) const {
            return (x == m.x && y == m.y);
        }

        bool operator!=(dim2 m) const {
            return !(x == m.x && y == m.y);
        }
    };

}
#endif
