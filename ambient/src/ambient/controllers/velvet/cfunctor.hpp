namespace ambient { namespace controllers { namespace velvet {

    constexpr int lim_expr(int n){ return n < AMBIENT_NUM_PROCS ? n : AMBIENT_NUM_PROCS; }
    constexpr int and_expr(bool m, bool n){ return (m && n); }

    template <int First, int Last>
    struct static_for {
         template <typename Fn>
         bool operator()(Fn const& fn) const {
             return First < Last ? and_expr(static_for<lim_expr(First+1), Last>()(fn), fn(First)) : true;
         }
    };

    inline void cfunctor::queue(cfunctor* d){
        this->deps.push_back(d);
    }

    // {{{ revision get/set

    template<int N>
    inline void get<revision, N>::spawn(revision& r){
        if(r.transfer == NULL){
            r.generator = r.transfer = new get(r);
            ambient::controller.queue((get*)r.transfer);
        }
    }

    template<int N>
    inline get<revision, N>::get(revision& r){
        /*if(ambient::rank()){
            if(r.valid()) printf("REPEATED RECV OF SZ %d\n", (int)r.extent);
            else printf("RECV OF SZ %d\n", (int)r.extent);
        }*/
        if(r.region != PERSIST_REGION)
            r.region = BULK_REGION;
        ambient::controller.alloc(r);
        handle = ambient::channel.get(&r);
        target = &r; 
    }

    template<int N>
    inline bool get<revision, N>::ready(){
        return ambient::channel.test(handle);
    }

    template<int N>
    inline void get<revision, N>::invoke(){
        target->complete();
        if(target->region != PERSIST_REGION)
            target->transfer = NULL;
    }

    template<int N>
    inline set<revision>& set<revision, N>::spawn(revision& r){
        if(r.transfer == NULL)
            r.transfer = new (r.region == PERSIST_REGION 
                              ? ambient::pool.malloc<sizeof(set)>()
                              : ambient::bulk.malloc<sizeof(set)>())
                         set<revision>(r);
        return *(set<revision>*)r.transfer;
    }

    template<int N>
    inline set<revision, N>::set(revision& r) : target(&r) {
        states = new std::vector<bool>(AMBIENT_NUM_PROCS, false);
        handles = new std::vector<request*>();
    }

    template<int N>
    template<int NE>
    inline set<revision, N>::set(set<revision, NE>* s) 
    : evaluated(false), target(s->target),
      states(s->states), handles(s->handles)
    {
    }

    template<int N>
    inline void set<revision, N>::operator >> (int p){
        if(!(*states)[p]){
            (*states)[p] = true;
            handles->push_back((request*)p);
            new (this) set<revision, lim_expr(N+1)>(this);
            if(!N){
                target->use();
                if(target->generator != NULL) ((cfunctor*)target->generator)->queue(this);
                else ambient::controller.queue(this);
            }
        }
    }

    template<int N>
    inline bool set<revision, N>::ready(){
        if(target->generator != NULL) return false;
        if(!evaluated){
            evaluated = true;
            for(int i = 0; i < N; ++i)
                (*handles)[i] = ambient::channel.set(target, (size_t)(*handles)[i]); 
            //if(N == (AMBIENT_NUM_PROCS-1)) printf("IT SHOULD BE A BROADCAST!\n");
        }
        return static_for<0, lim_expr(N)>()([&](int i){ 
            return ambient::channel.test((*handles)[i]);  // do we need to remove completed reqs?
        });
    }

    template<int N>
    inline void set<revision, N>::invoke(){
        target->release(); 
        if(target->region == PERSIST_REGION){
            this->handles->clear();
            new (this) set<revision>(this);
        }else{
            target->transfer = NULL;
            delete this->handles;
            delete this->states;
        }
    }

    // }}}

    // {{{ revision broadcast (get)/set

    inline void set<revision, AMBIENT_NUM_PROCS>::spawn(revision& r){
        if(r.transfer == NULL)
            r.transfer = new (r.region == PERSIST_REGION 
                              ? ambient::pool.malloc<sizeof(set)>()
                              : ambient::bulk.malloc<sizeof(set)>())
                         set<revision>(r);
        for(int i = 0; i < AMBIENT_NUM_PROCS; ++i) 
            (*(set<revision>*)r.transfer) >> i;
    }

    template<int NE>
    inline set<revision, AMBIENT_NUM_PROCS>::set(set<revision, NE>* s) 
    : evaluated(false), target(s->target),
      states(s->states), handles(s->handles)
    {
    }

    inline void set<revision, AMBIENT_NUM_PROCS>::operator >> (int p){ }

    inline bool set<revision, AMBIENT_NUM_PROCS>::ready(){
        if(target->generator != NULL) return false;
        if(!evaluated){ 
            evaluated = true; 
            for(int i = 0; i < AMBIENT_NUM_PROCS; ++i)
                (*handles)[i] = ambient::channel.set(target, (size_t)(*handles)[i]); 
        }
        return static_for<0, lim_expr(AMBIENT_NUM_PROCS)>()([&](int i){ 
            return ambient::channel.test((*handles)[i]);
        });
    }

    inline void set<revision, AMBIENT_NUM_PROCS>::invoke(){
        target->release(); 
        if(target->region == PERSIST_REGION){
            this->handles->clear();
            new (this) set<revision>(this);
        }else{
            target->transfer = NULL;
            delete this->handles;
            delete this->states;
        }
    }

    // }}}

    // {{{ transformable broadcast get/set

    inline void get<transformable, AMBIENT_NUM_PROCS>::spawn(transformable& v){
        ambient::controller.queue(new get(v));
    }

    inline get<transformable, AMBIENT_NUM_PROCS>::get(transformable& v){
        handle = ambient::channel.get(&v);
    }

    inline bool get<transformable, AMBIENT_NUM_PROCS>::ready(){
        return ambient::channel.test(handle);
    }

    inline void get<transformable, AMBIENT_NUM_PROCS>::invoke(){}

    inline set<transformable, AMBIENT_NUM_PROCS>& set<transformable, AMBIENT_NUM_PROCS>::spawn(transformable& v){
        set<transformable, AMBIENT_NUM_PROCS>* transfer = new set<transformable, AMBIENT_NUM_PROCS>(v);
        ((cfunctor*)v.generator)->queue(transfer);
        return *transfer;
    }

    inline set<transformable, AMBIENT_NUM_PROCS>::set(transformable& v) 
    : target(&v), evaluated(false) {
        handles = new std::vector<request*>(AMBIENT_NUM_PROCS);
    }

    inline bool set<transformable, AMBIENT_NUM_PROCS>::ready(){
        if(target->generator != NULL) return false;
        if(!evaluated){ 
            evaluated = true; 
            for(int i = 0; i < AMBIENT_NUM_PROCS; ++i)
                (*handles)[i] = ambient::channel.set(target, i); 
        }
        return static_for<0, lim_expr(AMBIENT_NUM_PROCS)>()([&](int i){ 
            return ambient::channel.test((*handles)[i]);
        });
    }

    inline void set<transformable, AMBIENT_NUM_PROCS>::invoke(){
        delete this->handles;
    }

    // }}}

} } }
