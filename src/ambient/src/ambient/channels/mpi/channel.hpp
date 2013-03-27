#include "ambient/channels/mpi/packets/auxiliary.hpp"
#include "ambient/utils/timings.hpp"

#define AMBIENT_MASTER_RANK 0
#define RESERVATION 2

namespace ambient { namespace controllers { namespace velvet {

//    void deviate(ambient::channels::mpi::packet&);

} } }

namespace ambient { namespace channels { namespace mpi {

    // {{{ channel

    inline channel::channel(){
        pthread_mutex_init(&this->mutex, NULL);
    }

    inline channel::~channel(){
        MPI_Finalize();
        pthread_mutex_destroy(&this->mutex);
    }

    inline void channel::init(){
        int level, zero = 0;
        MPI_Init_thread(&zero, NULL, MPI_THREAD_MULTIPLE, &level); 
        if(level != MPI_THREAD_MULTIPLE) printf("ERROR: Wrong threading level\n");
        this->world = new group(AMBIENT_MASTER_RANK, MPI_COMM_WORLD);
        this->volume = this->world->size;

//        this->add_handler( get_t<control_packet_t>() , controllers::velvet::deviate);
        //this->active = true;
        //pthread_create(&this->thread, NULL, &channel::stream, this);
    }

    inline void* channel::stream(void* instance){
        channel* c = static_cast<channel*>(instance);
        while(c->active) c->spin();
        return NULL;
    }

    inline void channel::spin(){
        pthread_mutex_lock(&this->mutex);
        std::list<pipe*>::const_iterator begin = this->qs.begin();
        std::list<pipe*>::const_iterator end = this->qs.end();
        pthread_mutex_unlock(&this->mutex);
        for(std::list<pipe*>::const_iterator p = begin; p != end; ++p)
        (*p)->spin();
    }

    inline void channel::add_handler(const packet_t& type, void(*callback)(packet&)){
        this->get_pipe(type, pipe::IN)->packet_delivered += callback;
    }

    inline channel::pipe* channel::get_pipe(const packet_t& type, channel::pipe::direction flow){
        for(std::list<pipe*>::const_iterator it = this->qs.begin(); it != qs.end(); ++it){
            if(&((*it)->type) == &(type) && flow == (*it)->flow){ return (*it); }
        }
        return this->add_pipe(type, flow); // create new if not found
    }

    inline channel::pipe* channel::add_pipe(const packet_t& type, channel::pipe::direction flow){
        pthread_mutex_lock(&this->mutex);
        this->qs.push_back(new pipe(type, flow));
        pthread_mutex_unlock(&this->mutex);
        return this->qs.back();
    }


    inline request::request(void* memory) 
    : memory(memory), mpi_request(MPI_REQUEST_NULL) 
    {
    };

    inline request* channel::get(transformable* v){
        request* q = new request(&v->v);
        MPI_Irecv(q->memory, 
                  (int)sizeof(transformable::numeric_union)/sizeof(double), 
                  MPI_DOUBLE, 
                  MPI_ANY_SOURCE, 
                  v->sid, 
                  MPI_COMM_WORLD, 
                  &(q->mpi_request));
        return q;
    }

    inline request* channel::set(transformable* v, int rank){
        if(rank == ambient::rank()) return NULL;
        request* q = new request(&v->v);
        MPI_Isend(q->memory, 
                  (int)sizeof(transformable::numeric_union)/sizeof(double), 
                  MPI_DOUBLE, 
                  rank, 
                  v->sid, 
                  MPI_COMM_WORLD, 
                  &(q->mpi_request));
        return q;
    }

    inline request* channel::get(revision* r){
        request* q = new request(r->header);
        //printf("Getting %d with tag %d\n", (int)r->extent, r->sid);
        MPI_Irecv(q->memory, 
                  (int)r->extent/sizeof(double), 
                  MPI_DOUBLE, 
                  MPI_ANY_SOURCE, 
                  r->sid, 
                  MPI_COMM_WORLD, 
                  &(q->mpi_request));
        return q;
    }

    inline request* channel::set(revision* r, int rank){
        request* q = new request(r->header);
        //printf("Sending %d with tag %d to %d\n", (int)r->extent, r->sid, rank);
        MPI_Isend(q->memory, 
                  (int)r->extent/sizeof(double), 
                  MPI_DOUBLE, 
                  rank, 
                  r->sid, 
                  MPI_COMM_WORLD, 
                  &(q->mpi_request));
        return q;
    }

    inline bool channel::test(request* q){
        if(q == NULL) return true; // for transformable
        int flag = 0;
        MPI_Test(&(q->mpi_request), &flag, MPI_STATUS_IGNORE);
        if(flag) return true;
        return false;
    }

    inline void channel::replicate(vbp& p){
        //pipe* queue = this->get_pipe(p->get_t(), pipe::OUT);
        //queue->send(queue->attach_request(&p));
    }

    inline void channel::broadcast(vbp& p, bool root){
        if(!root){ replicate(p); return; }
        for(int i = 0; i < world->size; i++){
            if(i == ambient::rank()) continue;
            p.dest = (double)i;
            replicate(p);
        }
    }

    inline void channel::emit(packet* p){
        packet* pk = static_cast<packet*>(p);
        if(pk->get<int>(A_DEST_FIELD) == ambient::rank()){
            this->get_pipe(pk->get_t(), pipe::IN)->packet_delivered(p);
        }else{
            pipe* queue = this->get_pipe(p->get_t(), pipe::OUT);
            queue->send(queue->attach_request(pk->data));
        }
    }

    inline void channel::ifetch(group* g, size_t sid, size_t x, size_t y){
        for(int i = 0; i < g->size; i++){
            this->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), // pack from auxiliary
                                                       g->get_proc(i), "P2P", 
                                                       "INFORM OWNER ABOUT REQUEST",
                                                       sid, sid, "GENERIC", // fixme
                                                       ambient::rank(), // forward target
                                                       x, y));
        }
    }

    // }}}

    // {{{ channel::pipe

    inline channel::pipe::request::request(void* memory) 
    : memory(memory), mpi_request(MPI_REQUEST_NULL) 
    {
    };

    inline channel::pipe::pipe(const packet_t& type, channel::pipe::direction flow) 
    : type(type), flow(flow), packet_delivered() 
    {
        pthread_mutex_init(&this->reqs_mutex, NULL);
        if(flow == IN) for(int i=0; i < RESERVATION; i++)
            this->recv(this->create_request());
    }

    inline channel::pipe::~pipe(){ 
        pthread_mutex_destroy(&this->reqs_mutex);
        /* cancelling requests here */
    }

    inline size_t channel::pipe::get_bound() const {
        return (size_t)this->type.displacements[A_BLOCK_P_DATA_FIELD];
    }

    inline channel::pipe::request* channel::pipe::create_request(){
    // used only in constructor, don't need mutex (as guarded by qs_mutex)
        this->reqs.push_back(new request(alloc_t(this->type)));
        return this->reqs.back();
    }

    inline channel::pipe::request* channel::pipe::attach_request(void* memory){
        pthread_mutex_lock(&this->reqs_mutex);
        this->reqs.push_back(new request(NULL));
        request* target = this->reqs.back();
        target->memory = memory;
        pthread_mutex_unlock(&this->reqs_mutex);
        return target;
    }

    inline void channel::pipe::send(request* r){
        packet* p = unpack(this->type, r->memory);
        if(p->get<char>(A_OP_T_FIELD) == 'P'){
            MPI_Isend(r->memory, 1, this->type.mpi_t, p->get<int>(A_DEST_FIELD), this->type.t_code, MPI_COMM_WORLD, &(r->mpi_request));
        }else if(p->get<char>(A_OP_T_FIELD) == 'B'){ // unrolling into a series
            for(size_t i=0; i < (size_t)ambient::channel.world->size; i++){
                packet it = packet(*p);
                it.set(A_OP_T_FIELD, "P2P");
                it.set(A_DEST_FIELD, i);
                this->send(this->attach_request(it.data));
            }
        }
        delete p;
    }

    inline void channel::pipe::recv(request* r){
        MPI_Irecv(r->memory, 1, this->type.mpi_t, MPI_ANY_SOURCE, this->type.t_code, MPI_COMM_WORLD, &(r->mpi_request));
    }

    inline void channel::pipe::renew(request* r){
        r->memory = alloc_t(this->type);
        this->recv(r);
    }

    inline void channel::pipe::spin(){ // ex: stream thread
        int flag = 0;
        std::list<request*>::iterator it;
        pthread_mutex_lock(&this->reqs_mutex);
        std::list<request*>::iterator begin = this->reqs.begin();
        std::list<request*>::iterator end = this->reqs.end();
        pthread_mutex_unlock(&this->reqs_mutex);

        for(it = begin; it != end; ++it){
            MPI_Test(&((*it)->mpi_request), &flag, MPI_STATUS_IGNORE);
            if(flag){
                packet* p = unpack(this->type, (*it)->memory); 
                this->packet_delivered(*p); // delegate
                // checked_free(p->data); // freeing memory
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

} } }
