#ifndef AMBIENT_UTILS_MEMORY
#define AMBIENT_UTILS_MEMORY
#include "ambient/utils/singleton.hpp"
#include "boost/pool/singleton_pool.hpp"

#define FUTURE_SIZE 16
#define BULK_CHUNK 52428800
#define CPU __cilkrts_get_worker_number()

namespace ambient { namespace memory {

    class bulk : public singleton< bulk > 
    {
    private:
        template<size_t R>
        struct pool {
            pool(){
                this->memory = (void*)std::malloc(R);
            }
            void use(){
                this->iterator = (char*)this->memory;
            }
            template<size_t S>
            void* malloc(){
                return ((char(*&)[S])this->iterator)++;
            }
            char* iterator;
            void* memory;
        }; 
        
        typedef pool<BULK_CHUNK> chunk;
    public:
        bulk(){
            this->arity = __cilkrts_get_nworkers();
            this->pools = new std::vector<chunk*>[arity];
            this->iterators = (size_t*)calloc(arity, sizeof(size_t));
            for(int i = 0; i < this->arity; i++){
                this->pools[i].push_back(new chunk());
                this->pools[i][0]->use();
            }
        }
        template<size_t S>
        void* malloc(){
            size_t& i = this->iterators[CPU];
            std::vector<chunk*>& v = this->pools[CPU];

            if(((size_t)v[i]->iterator - (size_t)v[i]->memory + S) >= (size_t)BULK_CHUNK){
                if(++i == v.size()) v.push_back(new chunk());
                v[i]->use();
            }
            return v[i]->malloc<S>();
        }
        void drop(){
            memset(this->iterators, 0, sizeof(size_t)*arity);
            for(int i = 0; i < this->arity; i++) this->pools[CPU][0]->use();
        }
    private:
        std::vector<chunk*>* pools;
        size_t* iterators;
        size_t  arity;
    };


    class pool : public singleton< pool >
    {
    private:
        template<size_t S, size_t R>
        class heap {
            typedef char(bulk)[S];
        public:
            heap() : capacity(R) {
                this->memory.reserve(2);
                this->handles.reserve(R);
                this->realloc();
            }
           ~heap(){
                for(int i = 0; i < this->memory.size(); i++)
                    free(this->memory[i]);
            }
            void realloc(){
                bulk* memory = (bulk*)std::malloc(sizeof(bulk)*this->capacity);
                for(int i = 0; i < capacity; ++i) this->handles.push_back(&memory[i]);
                this->memory.push_back(memory);
                this->capacity *= 2;
            }
            void* malloc(){
                if(this->handles.empty()) this->realloc();
                void* entry = this->handles.back();
                this->handles.pop_back();
                return entry;
            }
            void free(void* ptr){
                this->handles.push_back(ptr);
            }
        private:
            size_t capacity;
            std::vector<void*> handles;
            std::vector<void*> memory;
        };
        struct empty { };
        typedef heap<3000,  1000>  heap_3000;  heap_3000* p3000;
        typedef heap<4000,  1000>  heap_4000;  heap_4000* p4000;
        typedef heap<5000,  1000>  heap_5000;  heap_5000* p5000;
        typedef heap<6000,  500>   heap_6000;  heap_6000* p6000;
        typedef heap<7000,  300>   heap_7000;  heap_7000* p7000;
        typedef heap<8000,  30 >   heap_8000;  heap_8000* p8000;
        typedef heap<9000,  100>   heap_9000;  heap_9000* p9000;
        typedef heap<10000, 300>   heap_10000; heap_10000* p10000;
        typedef heap<11000, 500>   heap_11000; heap_11000* p11000;
        typedef heap<12000, 100>   heap_12000; heap_12000* p12000;
        typedef heap<13000, 60 >   heap_13000; heap_13000* p13000;
        typedef heap<14000, 100>   heap_14000; heap_14000* p14000;
        typedef heap<15000, 60 >   heap_15000; heap_15000* p15000;
        typedef heap<16000, 30 >   heap_16000; heap_16000* p16000;
        typedef heap<17000, 30 >   heap_17000; heap_17000* p17000;
        typedef heap<18000, 100>   heap_18000; heap_18000* p18000;
        typedef heap<19000, 30 >   heap_19000; heap_19000* p19000;
        typedef heap<20000, 30 >   heap_20000; heap_20000* p20000;
        typedef heap<21000, 100>   heap_21000; heap_21000* p21000;
        typedef heap<22000, 300>   heap_22000; heap_22000* p22000;
        typedef heap<23000, 30 >   heap_23000; heap_23000* p23000;
        typedef heap<24000, 200>   heap_24000; heap_24000* p24000;
        typedef heap<25000, 300>   heap_25000; heap_25000* p25000;
        typedef heap<26000, 1600>  heap_26000; heap_26000* p26000;
        typedef heap<27000, 30 >   heap_27000; heap_27000* p27000;
        typedef heap<28000, 100>   heap_28000; heap_28000* p28000;
        typedef heap<29000, 30 >   heap_29000; heap_29000* p29000;
        typedef heap<30000, 30 >   heap_30000; heap_30000* p30000;
        typedef heap<31000, 30 >   heap_31000; heap_31000* p31000;
        typedef heap<32000, 60 >   heap_32000; heap_32000* p32000;
        typedef heap<33000, 100>   heap_33000; heap_33000* p33000;
        typedef heap<34000, 100>   heap_34000; heap_34000* p34000;
        typedef heap<35000, 30 >   heap_35000; heap_35000* p35000;
        typedef heap<36000, 200>   heap_36000; heap_36000* p36000;
        typedef heap<37000, 30 >   heap_37000; heap_37000* p37000;
        typedef heap<38000, 30 >   heap_38000; heap_38000* p38000;
        typedef heap<39000, 30 >   heap_39000; heap_39000* p39000;
        typedef heap<40000, 100>   heap_40000; heap_40000* p40000;
        typedef heap<41000, 30 >   heap_41000; heap_41000* p41000;
        typedef heap<42000, 30 >   heap_42000; heap_42000* p42000;
        typedef heap<43000, 30 >   heap_43000; heap_43000* p43000;
        typedef heap<44000, 30 >   heap_44000; heap_44000* p44000;
        typedef heap<45000, 30 >   heap_45000; heap_45000* p45000;
        typedef heap<46000, 100>   heap_46000; heap_46000* p46000;
        typedef heap<47000, 100>   heap_47000; heap_47000* p47000;
        typedef heap<48000, 30 >   heap_48000; heap_48000* p48000;
        typedef heap<49000, 90 >   heap_49000; heap_49000* p49000;
        typedef heap<50000, 30 >   heap_50000; heap_50000* p50000;
        typedef heap<51000, 30 >   heap_51000; heap_51000* p51000;
        typedef heap<52000, 30 >   heap_52000; heap_52000* p52000;
        typedef heap<53000, 30 >   heap_53000; heap_53000* p53000;
        typedef heap<54000, 30 >   heap_54000; heap_54000* p54000;
        typedef heap<55000, 30 >   heap_55000; heap_55000* p55000;
        typedef heap<56000, 30 >   heap_56000; heap_56000* p56000;
        typedef heap<57000, 30 >   heap_57000; heap_57000* p57000;
        typedef heap<58000, 30 >   heap_58000; heap_58000* p58000;
        typedef heap<59000, 30 >   heap_59000; heap_59000* p59000;
        typedef heap<60000, 30 >   heap_60000; heap_60000* p60000;
        typedef heap<61000, 30 >   heap_61000; heap_61000* p61000;
        typedef heap<62000, 30 >   heap_62000; heap_62000* p62000;
        typedef heap<63000, 30 >   heap_63000; heap_63000* p63000;
        typedef heap<64000, 30 >   heap_64000; heap_64000* p64000;
        typedef heap<65000, 30 >   heap_65000; heap_65000* p65000;
        typedef heap<66000, 30 >   heap_66000; heap_66000* p66000;

        typedef heap<2097152, 30 >   heap_IB; heap_IB* pIB;
    public:
        pool(){
            this->arity = __cilkrts_get_nworkers();
            this->pIB   = new heap_IB  [arity];
        }
       ~pool(){
            delete[] this->pIB;
        }
        void* malloc(size_t sz){
            int ts = CPU;
            if(sz == 2097152) return pIB[ts].malloc();

            switch((int)((sz-1) / 1000)){
                    default: return std::malloc(sz); 
            }
        }
        template<typename T>
        static void* malloc(){
             return boost::singleton_pool<empty, sizeof(T)>::malloc();
        }
        template<size_t S>
        static void* malloc(){
             return boost::singleton_pool<empty, S>::malloc();
        }
        void free(void* ptr, size_t sz){
            if(ptr == NULL) return;
            int ts = (int)(drand48()*(double)this->arity);
            if(sz == 2097152) return pIB[ts].free(ptr);

            switch((int)((sz-1) / 1000)){
                    default: return std::free(ptr);
            }
        }
        template<typename T>
        static void free(void* ptr){
             boost::singleton_pool<empty, sizeof(T)>::free(ptr);
        }
        template<size_t S>
        static void free(void* ptr){
             boost::singleton_pool<empty, S>::free(ptr);
        }
    private:
        int arity;
    };

} }

namespace ambient {
    extern memory::bulk& bulk;
    extern memory::pool& pool;
}

#undef CPU
#endif
