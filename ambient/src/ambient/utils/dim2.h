#ifndef AMBIENT_UTILS_DIM2
#define AMBIENT_UTILS_DIM2

namespace ambient{

    class dim2 {
    public:
        size_t x, y;
        dim2() 
        : x(0), y(0) 
        {
        }

        dim2(size_t x, size_t y) 
        : x(x), y(y) 
        {
        }

        size_t max() const {
            return std::max(this->x, this->y);
        }

        size_t min() const {
            return std::min(this->x, this->y);
        }

        size_t square() const {
            return this->x * this->y;
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

        bool operator!=(int value) const {
            return !(x == value && y == value);
        }

        bool operator==(dim2 m) const {
            return (x == m.x && y == m.y);
        }

        bool operator!=(dim2 m) const {
            return !(x == m.x && y == m.y);
        }

        bool operator<(dim2 m) const {
            return (x < m.x && y < m.y);
        }

        bool operator<(size_t m) const {
            return (x < m && y < m);
        }

        bool operator>(dim2 m) const {
            return (x > m.x && y > m.y);
        }

        bool operator>(size_t m) const {
            return (x > m && y > m);
        }
    };

}

#endif
