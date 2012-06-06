#ifndef AMBIENT_UTILS_DIM2
#define AMBIENT_UTILS_DIM2

namespace ambient{

    class dim2 {
    public:
        size_t x, y;
        inline dim2() 
        : x(0), y(0) 
        {
        }

        inline dim2(size_t x, size_t y) 
        : x(x), y(y) 
        {
        }

        inline size_t max(){
            return std::max(this->x, this->y);
        }

        inline size_t min(){
            return std::min(this->x, this->y);
        }

        inline size_t square(){
            return this->x * this->y;
        }

        inline dim2& operator=(int value){
            this->x = this->y = value;
            return *this;
        }

        inline dim2& operator*=(const dim2 & b){
            this->x *= b.x;
            this->y *= b.y;
            return *this;
        }

        inline size_t operator*(const dim2 & b) const { // multiplication of all components
            return this->x * b.x *
                   this->y * b.y ;
        }

        inline bool operator==(int value) const {
            return (x == value && y == value);
        }

        inline bool operator!=(int value) const {
            return !(x == value && y == value);
        }

        inline bool operator==(dim2 m) const {
            return (x == m.x && y == m.y);
        }

        inline bool operator!=(dim2 m) const {
            return !(x == m.x && y == m.y);
        }

        inline bool operator<(dim2 m) const {
            return (x < m.x && y < m.y);
        }

        inline bool operator<(size_t m) const {
            return (x < m && y < m);
        }

        inline bool operator>(dim2 m) const {
            return (x > m.x && y > m.y);
        }

        inline bool operator>(size_t m) const {
            return (x > m && y > m);
        }
    };

}

#endif
