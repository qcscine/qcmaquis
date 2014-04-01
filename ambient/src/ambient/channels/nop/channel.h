/*
 * Ambient, License - Version 1.0 - May 3rd, 2012
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

#ifndef AMBIENT_CHANNELS_NOP_CHANNEL
#define AMBIENT_CHANNELS_NOP_CHANNEL

#define AMBIENT_CHANNEL_NAME nop

namespace ambient { namespace channels { namespace nop {

    class multirank {
    public:
        int operator()() const { return 0; }
        int left_neighbor() const { return 0; }
        int right_neighbor() const { return 0; }
    };

    template<class T> struct collective {
        bool test(){ return true; }
        void operator += (int rank){}
        bool involved(){ return true; }
    };

    class channel {
    public:
        typedef typename ambient::models::ssm::revision block_type;
        typedef typename ambient::models::ssm::transformable scalar_type;
        template<class T> using collective_type = collective<T>;
        size_t dim() const { return 1; }
        static void barrier(){}
        collective<block_type>* get(block_type& r){ return NULL; }
        collective<block_type>* set(block_type& r){ return NULL; }
        collective<scalar_type>* bcast(scalar_type& v, int root){ return NULL; }
        collective<scalar_type>* bcast(scalar_type& v){ return NULL; }
        multirank rank;
    };

} } }

#endif
