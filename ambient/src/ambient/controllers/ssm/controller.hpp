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

#include "ambient/utils/io.hpp"
#include "ambient/utils/timings.hpp"
#include "ambient/utils/overseer.hpp"
#include "ambient/utils/service.hpp"
#include "ambient/utils/mem.h"

namespace ambient { namespace controllers { namespace ssm {

    inline controller::~controller(){ 
        if(!chains->empty()) printf("Ambient:: exiting with operations still in queue!\n");
        this->clear();
    }

    inline controller::controller(){
        this->sid = 13;
        this->stack_m.reserve(AMBIENT_STACK_RESERVE);
        this->stack_s.reserve(AMBIENT_STACK_RESERVE);
        this->chains = &this->stack_m;
        this->mirror = &this->stack_s;
    }

    inline void controller::init(int db){
        this->db = get_num_procs() > db ? db : 0;
        this->context_base = new ambient::scope<base>();
        this->context = this->context_base;
        this->serial = (get_num_procs() == 1) ? true : false;
    }

    inline bool controller::tunable(){
        if(serial) return false;
        return context->tunable();
    }

    inline void controller::schedule(){
        const_cast<scope*>(context)->toss();
    }

    inline void controller::intend_read(history* o){
        revision* r = o->back(); if(r == NULL || model.common(r)) return;
        int candidate = model.remote(r) ? r->owner : (int)ambient::rank();
        context->score(candidate, r->spec.extent);
    }

    inline void controller::intend_write(history* o){
        revision* r = o->back(); if(r == NULL || model.common(r)) return;
        int candidate = model.remote(r) ? r->owner : (int)ambient::rank();
        context->select(candidate);
    }

    inline bool controller::scoped() const {
        return (context != context_base);
    }

    inline void controller::set_context(const scope* s){
        this->context = s; // no nesting
    }

    inline void controller::pop_context(){
        this->context = this->context_base;
    }

    inline bool controller::remote(){
        return (this->context->state == ambient::remote);
    }

    inline bool controller::local(){
        return (this->context->state == ambient::local);
    }

    inline bool controller::common(){
        return (this->context->state == ambient::common);
    }

    inline int controller::which(){
        return this->context->rank;
    }

    inline void controller::flush(){
        typedef typename std::vector<functor*>::const_iterator veci;
        AMBIENT_SMP_ENABLE
        while(!chains->empty()){
            for(veci i = chains->begin(); i != chains->end(); ++i){
                if((*i)->ready()){
                    functor* task = *i;
                    AMBIENT_THREAD task->invoke();
                    int size = task->deps.size();
                    for(int n = 0; n < size; n++) task->deps[n]->ready();
                    mirror->insert(mirror->end(), task->deps.begin(), task->deps.end());
                }else mirror->push_back(*i);
            }
            chains->clear();
            std::swap(chains,mirror);
        }
        AMBIENT_SMP_DISABLE
        model.clock++;
        fence();
    }

    inline bool controller::empty(){
        return this->chains->empty();
    }

    inline void controller::clear(){
        this->garbage.clear();
    }

    inline bool controller::queue(functor* f){
        this->chains->push_back(f);
        return true;
    }

    inline bool controller::update(revision& r){
        if(r.assist.first != model.clock){
            r.assist.first = model.clock;
            return true;
        }
        return false;
    }

    inline void controller::sync(revision* r){
        if(serial) return;
        if(model.common(r)) return;
        if(model.feeds(r)) set<revision>::spawn(*r);
        else get<revision>::spawn(*r);
    }

    inline void controller::lsync(revision* r){
        if(model.common(r)) return;
        if(!model.feeds(r)) get<revision>::spawn(*r);
    }

    inline void controller::rsync(revision* r){
        if(model.common(r)) return;
        if(r->owner != which()){
            if(model.feeds(r)) set<revision>::spawn(*r);
            else get<revision>::spawn(*r); // assist
        }
    }

    inline void controller::lsync(transformable* v){
        if(serial) return;
        set<transformable>::spawn(*v);
    }

    inline void controller::rsync(transformable* v){
        get<transformable>::spawn(*v);
    }

    template<typename T> void controller::collect(T* o){
        this->garbage.push_back(o);
    }

    inline void controller::squeeze(revision* r) const {
        if(r->valid() && !r->referenced() && r->locked_once()){
            if(r->spec.region == ambient::rstandard){
                ambient::pool::free(r->data, r->spec);
                r->spec.region = ambient::rdelegated;
            }else if(r->spec.region == ambient::rbulked){
                ambient::memory::data_bulk::reuse(r->data);
                r->spec.region = ambient::rdelegated;
            }
        }
    }

    inline void controller::touch(const history* o){
        model.touch(o);
    }

    inline void controller::use_revision(history* o){
        model.use_revision(o);
    }

    template<ambient::locality L, typename G>
    void controller::add_revision(history* o, G g){
        model.add_revision<L>(o, g);
    }

    inline int controller::get_rank() const {
        return channel.rank();
    }

    inline int controller::get_shared_rank() const {
        return get_num_procs();
    }

    inline int controller::get_dedicated_rank() const {
        return (get_num_procs()-db);
    }
        
    inline bool controller::verbose() const {
        return (get_rank() == 0);
    }

    inline void controller::fence() const {
        channel.barrier();
    }

    inline void controller::meminfo() const {
        size_t current_size = getCurrentRSS();
        size_t peak_size = getPeakRSS();
        for(int i = 0; i < get_num_procs(); i++){
            if(get_rank() == i){
                std::cout << "R" << i << ": current size: " << current_size << " " << current_size/1024/1024/1024 << "G.\n";
                std::cout << "R" << i << ": peak size: " << peak_size << " " << peak_size/1024/1024/1024 << "G.\n";
            }
            fence();
        }
    }

    inline int controller::get_num_procs() const {
        return channel.dim();
    }

    inline int controller::get_num_db_procs() const {
        return this->db;
    }

    inline int controller::get_num_workers() const {
        return (get_num_procs()-db);
    }

    inline int controller::generate_sid(){
        ++this->sid %= AMBIENT_MAX_SID; 
        return this->sid;
    }

    inline int controller::get_sid() const {
        return this->sid;
    }

    inline typename controller::channel_type & controller::get_channel(){
        return channel;
    }

} } }
