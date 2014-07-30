/*
 * Ambient Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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
