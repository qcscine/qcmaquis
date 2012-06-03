#ifndef AMBIENT_CHANNELS_MPI_GROUPS_MULTIRANK
#define AMBIENT_CHANNELS_MPI_GROUPS_MULTIRANK

#define UNDEFINED_RANK MPI_UNDEFINED

namespace ambient { namespace channels { namespace mpi {

    class multirank : public singleton< multirank >
    {
    public:
        inline multirank(){
        }
    private: 
        std::map<std::string,int> map;
    public:
        inline int operator()(const group* grp) const { return grp->rank; }
        inline int operator()(const char* name = "ambient") const {
            group* grp = group_map(name);
            if(grp != NULL) return grp->rank;
            return MPI_UNDEFINED;
        }
    };

} } }

namespace ambient {
    extern channels::mpi::multirank& rank;
}

#endif
