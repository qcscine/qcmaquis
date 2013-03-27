#define UNDEFINED_RANK MPI_UNDEFINED

namespace ambient { namespace channels { namespace mpi {

    inline int multirank::operator()() const {
        return ambient::channel.world->rank;
    }

    inline int multirank::operator()(const group* grp) const { 
        return grp->rank; 
    }

    inline int multirank::translate(int rank, const group* source) const {
        if(source->depth == 0) return rank;
        return cast_to_parent(rank, source, ambient::channel.world);
    }

    inline int multirank::translate(int rank, const group* source, const group* target) const {
        if(source->depth == target->depth) return rank;
        else if(target->depth < source->depth) return cast_to_parent(rank, source, target);
        else return cast_to_child(rank, source, target);
    }

    // query rank "i" inside a group (not ranks[i])
    inline int multirank::cast_to_parent(int rank, const group* source, const group* target) const {
        for(const group* i = source; i != target; i = i->parent){
            assert(rank < i->size);
            rank = i->ranks[rank];
        }
        return rank;
    }
 
    inline int multirank::cast_to_child(int rank, const group* source, const group* target) const {
        if(target == source) return rank;
        rank = cast_to_child(rank, source, target->parent);
        for(int i = 0; i < target->size; ++i)
            if(target->ranks[i] == rank) return i;
        return UNDEFINED_RANK;
    }

    inline bool multirank::belongs(const group* target) const {
        return (target->rank != UNDEFINED_RANK);
    }

    inline bool multirank::masters(const group* target) const {
        return (target->rank == target->master);
    }

    inline void multirank::mute(){
        this->verbose = false;
    }

} } }
