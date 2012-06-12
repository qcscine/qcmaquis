#include "ambient/channels/mpi/packets/auxiliary.hpp"
#include "ambient/utils/timings.hpp"

#define AMBIENT_MASTER_RANK 0
#define RESERVATION 2

namespace ambient { namespace controllers { namespace velvet {

    void forward_block(ambient::channels::mpi::packet&);
    void accept_block(ambient::channels::mpi::packet&);

} } }

namespace ambient { namespace channels { namespace mpi {

    // {{{ channel

    inline channel::channel(){ 
        pthread_mutex_init(&this->mutex, NULL);
    }

    inline channel::~channel(){
        if(this->active) this->finalize();
        pthread_mutex_destroy(&this->mutex);
    }

    inline group* channel::world(){
        return this->ambient;
    }

    inline void channel::init(){
        int threading_level;
        MPI_Init_thread(0, NULL, MPI_THREAD_MULTIPLE, &threading_level);
        if(threading_level != MPI_THREAD_MULTIPLE) printf("Wrong value of threading_level!\n");
        assert(threading_level == MPI_THREAD_MULTIPLE);
        this->ambient = new group("ambient", AMBIENT_MASTER_RANK, MPI_COMM_WORLD);

        this->add_handler( get_t<layout_packet_t>() , controllers::velvet::forward_block );

        this->active = true;
        //pthread_create(&this->thread, NULL, &channel::stream, this);
    }

    inline std::pair<size_t*,size_t> channel::id(){
        return this->ambient->id;
    }

    inline packet_t& channel::get_block_packet_type(size_t len){
        printf("WARNING: Getting packet type for some reason!\n");
        static packet_t* types[65538] = { NULL };
        size_t idx = len / 8; // for double
        if(types[idx] == NULL){
            types[idx] = new block_packet_t(len);
            //pt->commit(); // serial
            //this->add_handler(*pt, controllers::accept_block);
        }
        return *types[idx];
    }

    inline void channel::finalize(){
        this->active = false;
        //pthread_join(this->thread, NULL);
        MPI_Finalize();
    }

    inline void* channel::stream(void* instance){ // pthread bootstrapper
        channel* c = static_cast<channel*>(instance);
        while(c->active) c->spin();
        return NULL;
    }

    inline size_t channel::get_volume() const {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        return (size_t)size;
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

    inline void channel::emit(packet* p){
        packet* pk = static_cast<packet*>(p);
        if(pk->get<int>(A_DEST_FIELD) == ambient::rank()){
            this->get_pipe(pk->get_t(), pipe::IN)->packet_delivered(p);
        }else{
            pipe* queue = this->get_pipe(p->get_t(), pipe::OUT);
            queue->send(queue->attach_request(pk->data));
        }
    }

    inline void channel::spin(){ // ex: stream thread
        pthread_mutex_lock(&this->mutex);
        std::list<pipe*>::const_iterator begin = this->qs.begin();
        std::list<pipe*>::const_iterator end = this->qs.end();
        pthread_mutex_unlock(&this->mutex);
        for(std::list<pipe*>::const_iterator p = begin; p != end; ++p){
            (*p)->spin();
        }
    }

    inline void channel::ifetch(group* g, size_t sid, size_t x, size_t y){
        for(int i = 0; i < g->get_size(); i++){
            this->emit(pack<layout_packet_t>(alloc_t<layout_packet_t>(), // pack from auxiliary
                                                       g->get_member(i), "P2P", 
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
        use(r->memory); // using the memory until completion
        if(p->get<char>(A_OP_T_FIELD) == 'P'){
            MPI_Isend(r->memory, 1, this->type.mpi_t, p->get<int>(A_DEST_FIELD), this->type.t_code, MPI_COMM_WORLD, &(r->mpi_request));
        }else if(p->get<char>(A_OP_T_FIELD) == 'B'){ // unrolling into a series
            for(size_t i=0; i < ambient::channel.get_volume(); i++){
                packet it = packet(*p);
                it.set(A_OP_T_FIELD, "P2P");
                it.set(A_DEST_FIELD, i);
                this->send(this->attach_request(it.data));
            }
        }
        delete p;
    }

    inline void channel::pipe::recv(request* r){
        use(r->memory);
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

} } }
