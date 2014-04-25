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

namespace ambient { namespace models { namespace ssm {

    inline revision::revision(size_t extent, void* g, ambient::locality l, rank_t owner)
    : spec(extent), generator(g), state(l), 
      data(NULL), users(0), owner(owner)
    {
    }

    inline void revision::embed(void* ptr){
        data = ptr;
    }

    inline void revision::reuse(revision& r){
        data = r.data;
        spec.reuse(r.spec);
    }

    inline void revision::use(){
        ++users;
    }

    inline void revision::release(){
        --users;
    }

    inline void revision::complete(){
        generator = NULL;
    }

    inline void revision::invalidate(){
        data = NULL;
    }

    inline bool revision::locked() const {
        return (users != 0);
    }

    inline bool revision::locked_once() const {
        return (users == 1);
    }

    inline bool revision::valid() const {
        return (data != NULL);
    }

    inline bool revision::referenced() const {
        return (spec.crefs != 0);
    }

} } }
