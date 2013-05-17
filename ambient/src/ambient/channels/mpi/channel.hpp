#define AMBIENT_MASTER_RANK 0
#define AMBIENT_MPI_THREADING MPI_THREAD_FUNNELED // MPI_THREAD_MULTIPLE

namespace ambient { namespace channels { namespace mpi {

    inline channel::~channel(){
        MPI_Finalize();
    }

    inline void channel::init(){
        int level, zero = 0;
        MPI_Init_thread(&zero, NULL, AMBIENT_MPI_THREADING, &level);
        if(level != AMBIENT_MPI_THREADING) printf("ERROR: Wrong threading level\n");
        this->world = new group(AMBIENT_MASTER_RANK, MPI_COMM_WORLD);
        this->volume = this->world->size;
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
        request* q = new request(r->data);
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
        if(rank == ambient::rank()) return NULL;
        request* q = new request(r->data);
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

} } }
