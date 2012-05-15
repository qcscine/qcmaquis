#include "ambient/ambient.h"
#include "ambient/channels/mpi_channel.h"
#include "ambient/channels/packets/auxiliary.hpp"
#include "ambient/controllers/v_controller.h"
#include "ambient/utils/timings.h"

#define AMBIENT_MASTER_RANK 0
#define RESERVATION 2

namespace ambient { namespace channels {

    // {{{ mpi_channel

    mpi_channel::mpi_channel(){ 
        pthread_mutex_init(&this->mutex, NULL);
    }

    mpi_channel::~mpi_channel(){
        if(this->active) this->finalize();
        pthread_mutex_destroy(&this->mutex);
    }

    group* mpi_channel::world(){
        return this->ambient;
    }

    void mpi_channel::init(){
        int threading_level;
        MPI_Init_thread(0, NULL, MPI_THREAD_MULTIPLE, &threading_level);
        if(threading_level != MPI_THREAD_MULTIPLE) printf("Wrong value of threading_level!\n");
        assert(threading_level == MPI_THREAD_MULTIPLE);
        this->ambient = new group("ambient", AMBIENT_MASTER_RANK, MPI_COMM_WORLD);

        this->add_handler( get_t<layout_packet_t>() , controllers::forward_block );

        this->active = true;
        //pthread_create(&this->thread, NULL, &mpi_channel::stream, this);
    }

    std::pair<size_t*,size_t> mpi_channel::id(){
        return this->ambient->id;
    }

    channels::packet_t& mpi_channel::get_block_packet_type(size_t len){
        static std::map<size_t,channels::packet_t*> map;
        std::map<size_t,channels::packet_t*>::iterator value = map.find(len);
        channels::packet_t* pt;

        if(value != map.end()) 
            pt = (*value).second;
        else{
            pt = new channels::block_packet_t(len);
            //pt->commit(); // serial
            map.insert(std::pair<size_t,channels::packet_t*>(len,pt));
            //this->add_handler(*pt, controllers::accept_block);
        }
        return *pt;
    }

    void mpi_channel::finalize(){
        this->active = false;
        //pthread_join(this->thread, NULL);
        MPI_Finalize();
    }

    void* mpi_channel::stream(void* instance){ // pthread bootstrapper
        mpi_channel* channel = static_cast<channels::mpi_channel*>(instance);
        while(channel->active) channel->spin();
        return NULL;
    }

    size_t mpi_channel::get_volume() const {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        return (size_t)size;
    }

    void mpi_channel::add_handler(const packet_t& type, void(*callback)(channels::packet&)){
        this->get_pipe(type, pipe::IN)->packet_delivered += callback;
    }

    mpi_channel::pipe* mpi_channel::get_pipe(const packet_t& type, mpi_channel::pipe::direction flow){
        for(std::list<pipe*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
            if(&((*it)->type) == &(type) && flow == (*it)->flow){ return (*it); }
        }
        return this->add_pipe(type, flow); // create new if not found
    }

    mpi_channel::pipe* mpi_channel::add_pipe(const packet_t& type, mpi_channel::pipe::direction flow){
        pthread_mutex_lock(&this->mutex);
        this->qs.push_back(new pipe(type, flow));
        pthread_mutex_unlock(&this->mutex);
        return this->qs.back();
    }

    void mpi_channel::emit(channels::packet* p){
        channels::packet* pk = static_cast<channels::packet*>(p);
        if(pk->get<int>(A_DEST_FIELD) == ambient::rank()){
            this->get_pipe(pk->get_t(), pipe::IN)->packet_delivered(p);
        }else{
            pipe* queue = this->get_pipe(p->get_t(), pipe::OUT);
            queue->send(queue->attach_request(pk->data));
        }
    }

    void mpi_channel::spin(){ // ex: stream thread
        pthread_mutex_lock(&this->mutex);
        std::list<pipe*>::const_iterator begin = this->qs.begin();
        std::list<pipe*>::const_iterator end = this->qs.end();
        pthread_mutex_unlock(&this->mutex);
        for(std::list<pipe*>::const_iterator p = begin; p != end; ++p){
            (*p)->spin();
        }
    }

    void mpi_channel::ifetch(group* g, size_t gid, size_t sid, size_t i, size_t j){
        for(int i = 0; i < g->get_size(); i++){
            this->emit(channels::pack<layout_packet_t>(alloc_t<layout_packet_t>(), // pack from auxiliary
                                                       g->get_member(i), "P2P", 
                                                       "INFORM OWNER ABOUT REQUEST",
                                                       gid, sid, "GENERIC",
                                                       ambient::rank(), // forward target
                                                       i, j));
        }
    }

    // }}}

    // {{{ mpi_channel::pipe

    mpi_channel::pipe::request::request(void* memory) 
    : memory(memory), mpi_request(MPI_REQUEST_NULL) 
    {
    };

    mpi_channel::pipe::pipe(const packet_t& type, mpi_channel::pipe::direction flow) 
    : type(type), flow(flow), packet_delivered() 
    {
        pthread_mutex_init(&this->reqs_mutex, NULL);
        if(flow == IN) for(int i=0; i < RESERVATION; i++)
            this->recv(this->create_request());
    }

    mpi_channel::pipe::~pipe(){ 
        pthread_mutex_destroy(&this->reqs_mutex);
        /* cancelling requests here */
    }

    size_t mpi_channel::pipe::get_bound() const {
        return (size_t)this->type.displacements[A_BLOCK_P_DATA_FIELD];
    }

    mpi_channel::pipe::request* mpi_channel::pipe::create_request(){
    // used only in constructor, don't need mutex (as guarded by qs_mutex)
        this->reqs.push_back(new request(alloc_t(this->type)));
        return this->reqs.back();
    }

    mpi_channel::pipe::request* mpi_channel::pipe::attach_request(void* memory){
        pthread_mutex_lock(&this->reqs_mutex);
        this->reqs.push_back(new request(NULL));
        request* target = this->reqs.back();
        target->memory = memory;
        pthread_mutex_unlock(&this->reqs_mutex);
        return target;
    }

    void mpi_channel::pipe::send(request* r){
        channels::packet* p = unpack(this->type, r->memory);
        use(r->memory); // using the memory until completion
        if(p->get<char>(A_OP_T_FIELD) == 'P'){
            MPI_Isend(r->memory, 1, this->type.mpi_t, p->get<int>(A_DEST_FIELD), this->type.t_code, MPI_COMM_WORLD, &(r->mpi_request));
        }else if(p->get<char>(A_OP_T_FIELD) == 'B'){ // unrolling into a series
            for(size_t i=0; i < channel.get_volume(); i++){
                channels::packet it = channels::packet(*p);
                it.set(A_OP_T_FIELD, "P2P");
                it.set(A_DEST_FIELD, i);
                this->send(this->attach_request(it.data));
            }
        }
        delete p;
    }

    void mpi_channel::pipe::recv(request* r){
        use(r->memory);
        MPI_Irecv(r->memory, 1, this->type.mpi_t, MPI_ANY_SOURCE, this->type.t_code, MPI_COMM_WORLD, &(r->mpi_request));
    }

    void mpi_channel::pipe::renew(request* r){
        r->memory = alloc_t(this->type);
        this->recv(r);
    }

    void mpi_channel::pipe::spin(){ // ex: stream thread
        int flag = 0;
        std::list<request*>::iterator it;
        pthread_mutex_lock(&this->reqs_mutex);
        std::list<request*>::iterator begin = this->reqs.begin();
        std::list<request*>::iterator end = this->reqs.end();
        pthread_mutex_unlock(&this->reqs_mutex);

        for(it = begin; it != end; ++it){
            MPI_Test(&((*it)->mpi_request), &flag, MPI_STATUS_IGNORE);
            if(flag){
                channels::packet* p = unpack(this->type, (*it)->memory); 
                this->packet_delivered(*p); // delegate
                unuse(p->data); checked_free(p->data); // freeing memory
                delete p;

                if(this->flow == IN) this->renew(*it);
                else{
                    pthread_mutex_lock(&this->reqs_mutex);
                    it = this->reqs.erase(it);
                    pthread_mutex_unlock(&this->reqs_mutex);
                }
            }
        }
    }

    // }}}

} }
