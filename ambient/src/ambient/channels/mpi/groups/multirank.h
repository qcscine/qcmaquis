#ifndef AMBIENT_CHANNELS_MPI_GROUPS_MULTIRANK
#define AMBIENT_CHANNELS_MPI_GROUPS_MULTIRANK

namespace ambient { namespace channels { namespace mpi {

    class multirank : public singleton< multirank >
    {
    public:
        multirank() : verbose(true) {}
        int operator()() const;
        int operator()(const group* grp) const;
        int translate(int rank, const group* source) const; // default: world
        int translate(int rank, const group* source, const group* target) const;
        int cast_to_parent(int rank, const group* source, const group* target) const;
        int cast_to_child(int rank, const group* source, const group* target) const;
        bool belongs(const group* target) const;
        bool masters(const group* target) const;
        void mute();
        bool verbose;
    };

    //    Context misc functions:
    //
    //    int  get_master(){ return ambient::rank.translate(grp->master, grp); }
    //    bool involved()  { return ambient::rank.belongs(grp);                }
    //    bool is_master() { return ambient::rank.masters(grp);                }
    //    int  get_rank()  { return grp->rank;                                 }
    //    int  get_size()  { return grp->size;                                 }
    //

} } }

namespace ambient {
    extern channels::mpi::multirank& rank;
}

#endif
