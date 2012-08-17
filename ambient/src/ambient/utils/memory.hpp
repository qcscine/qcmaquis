#ifndef AMBIENT_UTILS_MEMORY
#define AMBIENT_UTILS_MEMORY
#include "ambient/utils/singleton.hpp"
#include "boost/pool/singleton_pool.hpp"

#define BULK_LENGTH 8388608*50
#define FUTURE_SIZE 16

namespace ambient { namespace utils {

    class static_memory {
    public:
        struct empty { };

        template<typename T>
        static void* malloc(){
             return boost::singleton_pool<empty, sizeof(T)>::malloc();
        }

        template<size_t S>
        static void* malloc(){
             return boost::singleton_pool<empty, S>::malloc();
        }

        template<typename T>
        static void free(void* ptr){
             boost::singleton_pool<empty, sizeof(T)>::free(ptr);
        }

        template<size_t S>
        static void free(void* ptr){
             boost::singleton_pool<empty, S>::free(ptr);
        }
    };

    class bulk_memory : public singleton< bulk_memory > 
    {
    public:
        bulk_memory(){
            this->arity = __cilkrts_get_nworkers();
            this->pools = (void**)malloc(sizeof(void*)*this->arity);
            for(int i = 0; i < this->arity; i++)
                this->pools[i] = malloc(BULK_LENGTH);
            this->iterators = (char**)malloc(sizeof(char*)*this->arity);
            for(int i = 0; i < this->arity; i++)
                this->iterators[i] = (char*)this->pools[i];
        }
       ~bulk_memory(){
            for(int i = 0; i < this->arity; i++)
                free(this->pools[i]);
            free(this->iterators);
            free(this->pools);
        }
        template<size_t S>
        void* get(){
            char*& iterator = this->iterators[__cilkrts_get_worker_number()];
            void* result = iterator;
            iterator += S; // 16*((size_t)(S/16)+1); // alignment variant
            return result;
        }
        void refresh(){
            for(int i = 0; i < this->arity; i++)
                this->iterators[i] = (char*)this->pools[i];
        }
    private:
        char** iterators;
        void** pools;
        int arity;
    };


    template<size_t S, size_t R>
    class pool {
    public:
        pool(){
            this->memory.reserve(2);
            this->handles.reserve(R);
            this->realloc();
        }
       ~pool(){
            for(int i = 0; i < this->memory.size(); i++)
                free(this->memory[i]);
        }
        void realloc(){
            this->memory.push_back(std::malloc(S*R));
            char* memory = (char*)this->memory.back();
            for(int i = 0; i < R; ++i)
                this->handles.push_back(memory + i*S);
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
        std::vector<void*> handles;
        std::vector<void*> memory;
    };

    class range_memory : public singleton< range_memory >
    {
    private:
        typedef pool<100,  15000> pool_100;  pool_100 * p100;
        typedef pool<200,  6000 > pool_200;  pool_200 * p200;
        typedef pool<300,  2500 > pool_300;  pool_300 * p300;
        typedef pool<400,  2500 > pool_400;  pool_400 * p400;
        typedef pool<500,  3000 > pool_500;  pool_500 * p500;
        typedef pool<600,  2500 > pool_600;  pool_600 * p600;
        typedef pool<700,  1500 > pool_700;  pool_700 * p700;
        typedef pool<800,  1000 > pool_800;  pool_800 * p800;
        typedef pool<900,  1000 > pool_900;  pool_900 * p900;
        typedef pool<1000, 1000 > pool_1000; pool_1000* p1000;
        typedef pool<1100, 1500 > pool_1100; pool_1100* p1100;
        typedef pool<1200, 1500 > pool_1200; pool_1200* p1200;
        typedef pool<1300, 1500 > pool_1300; pool_1300* p1300;
        typedef pool<1400, 500 >  pool_1400; pool_1400* p1400;
        typedef pool<1500, 2500 > pool_1500; pool_1500* p1500;
        typedef pool<1600, 10000> pool_1600; pool_1600* p1600;
        typedef pool<1700, 60 >   pool_1700; pool_1700* p1700;
        typedef pool<1800, 1500>  pool_1800; pool_1800* p1800;
        typedef pool<1900, 30 >   pool_1900; pool_1900* p1900;
        typedef pool<2000, 100 >  pool_2000; pool_2000* p2000;
        typedef pool<2800, 100 >  pool_2800; pool_2800* p2800;
        typedef pool<3200, 1000>  pool_3200; pool_3200* p3200;
        typedef pool<6400, 2500>  pool_6400; pool_6400* p6400;

        typedef pool<3000,  1000>  pool_3000;  pool_3000* p3000;
        typedef pool<4000,  1000>  pool_4000;  pool_4000* p4000;
        typedef pool<5000,  1000>  pool_5000;  pool_5000* p5000;
        typedef pool<6000,  500>   pool_6000;  pool_6000* p6000;
        typedef pool<7000,  300>   pool_7000;  pool_7000* p7000;
        typedef pool<8000,  30 >   pool_8000;  pool_8000* p8000;
        typedef pool<9000,  100>   pool_9000;  pool_9000* p9000;
        typedef pool<10000, 300>   pool_10000; pool_10000* p10000;
        typedef pool<11000, 500>   pool_11000; pool_11000* p11000;
        typedef pool<12000, 100>   pool_12000; pool_12000* p12000;
        typedef pool<13000, 60 >   pool_13000; pool_13000* p13000;
        typedef pool<14000, 100>   pool_14000; pool_14000* p14000;
        typedef pool<15000, 60 >   pool_15000; pool_15000* p15000;
        typedef pool<16000, 30 >   pool_16000; pool_16000* p16000;
        typedef pool<17000, 30 >   pool_17000; pool_17000* p17000;
        typedef pool<18000, 100>   pool_18000; pool_18000* p18000;
        typedef pool<19000, 30 >   pool_19000; pool_19000* p19000;
        typedef pool<20000, 30 >   pool_20000; pool_20000* p20000;
        typedef pool<21000, 100>   pool_21000; pool_21000* p21000;
        typedef pool<22000, 300>   pool_22000; pool_22000* p22000;
        typedef pool<23000, 30 >   pool_23000; pool_23000* p23000;
        typedef pool<24000, 200>   pool_24000; pool_24000* p24000;
        typedef pool<25000, 300>   pool_25000; pool_25000* p25000;
        typedef pool<26000, 1600>  pool_26000; pool_26000* p26000;
        typedef pool<27000, 30 >   pool_27000; pool_27000* p27000;
        typedef pool<28000, 100>   pool_28000; pool_28000* p28000;
        typedef pool<29000, 30 >   pool_29000; pool_29000* p29000;
        typedef pool<30000, 30 >   pool_30000; pool_30000* p30000;
        typedef pool<31000, 30 >   pool_31000; pool_31000* p31000;
        typedef pool<32000, 60 >   pool_32000; pool_32000* p32000;
        typedef pool<33000, 100>   pool_33000; pool_33000* p33000;
        typedef pool<34000, 100>   pool_34000; pool_34000* p34000;
        typedef pool<35000, 30 >   pool_35000; pool_35000* p35000;
        typedef pool<36000, 200>   pool_36000; pool_36000* p36000;
        typedef pool<37000, 30 >   pool_37000; pool_37000* p37000;
        typedef pool<38000, 30 >   pool_38000; pool_38000* p38000;
        typedef pool<39000, 30 >   pool_39000; pool_39000* p39000;
        typedef pool<40000, 100>   pool_40000; pool_40000* p40000;
        typedef pool<41000, 30 >   pool_41000; pool_41000* p41000;
        typedef pool<42000, 30 >   pool_42000; pool_42000* p42000;
        typedef pool<43000, 30 >   pool_43000; pool_43000* p43000;
        typedef pool<44000, 30 >   pool_44000; pool_44000* p44000;
        typedef pool<45000, 30 >   pool_45000; pool_45000* p45000;
        typedef pool<46000, 100>   pool_46000; pool_46000* p46000;
        typedef pool<47000, 100>   pool_47000; pool_47000* p47000;
        typedef pool<48000, 30 >   pool_48000; pool_48000* p48000;
        typedef pool<49000, 90 >   pool_49000; pool_49000* p49000;
        typedef pool<50000, 30 >   pool_50000; pool_50000* p50000;
        typedef pool<51000, 30 >   pool_51000; pool_51000* p51000;
        typedef pool<52000, 30 >   pool_52000; pool_52000* p52000;
        typedef pool<53000, 30 >   pool_53000; pool_53000* p53000;
        typedef pool<54000, 30 >   pool_54000; pool_54000* p54000;
        typedef pool<55000, 30 >   pool_55000; pool_55000* p55000;
        typedef pool<56000, 30 >   pool_56000; pool_56000* p56000;
        typedef pool<57000, 30 >   pool_57000; pool_57000* p57000;
        typedef pool<58000, 30 >   pool_58000; pool_58000* p58000;
        typedef pool<59000, 30 >   pool_59000; pool_59000* p59000;
        typedef pool<60000, 30 >   pool_60000; pool_60000* p60000;
        typedef pool<61000, 30 >   pool_61000; pool_61000* p61000;
        typedef pool<62000, 30 >   pool_62000; pool_62000* p62000;
        typedef pool<63000, 30 >   pool_63000; pool_63000* p63000;
        typedef pool<64000, 30 >   pool_64000; pool_64000* p64000;
        typedef pool<65000, 30 >   pool_65000; pool_65000* p65000;
        typedef pool<66000, 30 >   pool_66000; pool_66000* p66000;
    public:
        range_memory(){
            this->arity = __cilkrts_get_nworkers();
            this->p100  = new pool_100 [arity];
            this->p200  = new pool_200 [arity];
            this->p300  = new pool_300 [arity];
            this->p400  = new pool_400 [arity];
            this->p500  = new pool_500 [arity];
            this->p600  = new pool_600 [arity];
            this->p700  = new pool_700 [arity];
            this->p800  = new pool_800 [arity];
            this->p900  = new pool_900 [arity];
            this->p1000 = new pool_1000[arity];
            this->p1100 = new pool_1100[arity];
            this->p1200 = new pool_1200[arity];
            this->p1300 = new pool_1300[arity];
            this->p1400 = new pool_1400[arity];
            this->p1500 = new pool_1500[arity];
            this->p1600 = new pool_1600[arity];
            this->p1700 = new pool_1700[arity];
            this->p1800 = new pool_1800[arity];
            this->p1900 = new pool_1900[arity];
            this->p2000 = new pool_2000[arity];
            this->p2800 = new pool_2800[arity];
            this->p3200 = new pool_3200[arity];
            this->p6400 = new pool_6400[arity];

            this->p3000  = new pool_3000 [arity];
            this->p4000  = new pool_4000 [arity];
            this->p5000  = new pool_5000 [arity];
            this->p6000  = new pool_6000 [arity];
            this->p7000  = new pool_7000 [arity];
            this->p8000  = new pool_8000 [arity];
            this->p9000  = new pool_9000 [arity];
            this->p10000 = new pool_10000[arity];
            this->p11000 = new pool_11000[arity];
            this->p12000 = new pool_12000[arity];
            this->p13000 = new pool_13000[arity];
            this->p14000 = new pool_14000[arity];
            this->p15000 = new pool_15000[arity];
            this->p16000 = new pool_16000[arity];
            this->p17000 = new pool_17000[arity];
            this->p18000 = new pool_18000[arity];
            this->p19000 = new pool_19000[arity];
            this->p20000 = new pool_20000[arity];
            this->p21000 = new pool_21000[arity];
            this->p22000 = new pool_22000[arity];
            this->p23000 = new pool_23000[arity];
            this->p24000 = new pool_24000[arity];
            this->p25000 = new pool_25000[arity];
            this->p26000 = new pool_26000[arity];
            this->p27000 = new pool_27000[arity];
            this->p28000 = new pool_28000[arity];
            this->p29000 = new pool_29000[arity];
            this->p30000 = new pool_30000[arity];
            this->p31000 = new pool_31000[arity];
            this->p32000 = new pool_32000[arity];
            this->p33000 = new pool_33000[arity];
            this->p34000 = new pool_34000[arity];
            this->p35000 = new pool_35000[arity];
            this->p36000 = new pool_36000[arity];
            this->p37000 = new pool_37000[arity];
            this->p38000 = new pool_38000[arity];
            this->p39000 = new pool_39000[arity];
            this->p40000 = new pool_40000[arity];
            this->p41000 = new pool_41000[arity];
            this->p42000 = new pool_42000[arity];
            this->p43000 = new pool_43000[arity];
            this->p44000 = new pool_44000[arity];
            this->p45000 = new pool_45000[arity];
            this->p46000 = new pool_46000[arity];
            this->p47000 = new pool_47000[arity];
            this->p48000 = new pool_48000[arity];
            this->p49000 = new pool_49000[arity];
            this->p50000 = new pool_50000[arity];
            this->p51000 = new pool_51000[arity];
            this->p52000 = new pool_52000[arity];
            this->p53000 = new pool_53000[arity];
            this->p54000 = new pool_54000[arity];
            this->p55000 = new pool_55000[arity];
            this->p56000 = new pool_56000[arity];
            this->p57000 = new pool_57000[arity];
            this->p58000 = new pool_58000[arity];
            this->p59000 = new pool_59000[arity];
            this->p60000 = new pool_60000[arity];
            this->p61000 = new pool_61000[arity];
            this->p62000 = new pool_62000[arity];
            this->p63000 = new pool_63000[arity];
            this->p64000 = new pool_64000[arity];
            this->p65000 = new pool_65000[arity];
            this->p66000 = new pool_66000[arity];
        }
       ~range_memory(){
            delete[] this->p100;
            delete[] this->p200;
            delete[] this->p300;
            delete[] this->p400;
            delete[] this->p500;
            delete[] this->p600;
            delete[] this->p700;
            delete[] this->p800;
            delete[] this->p900;
            delete[] this->p1000;
            delete[] this->p1100;
            delete[] this->p1200;
            delete[] this->p1300;
            delete[] this->p1400;
            delete[] this->p1500;
            delete[] this->p1600;
            delete[] this->p1700;
            delete[] this->p1800;
            delete[] this->p1900;
            delete[] this->p2000;
            delete[] this->p2800;
            delete[] this->p3200;
            delete[] this->p6400;

            delete[] this->p3000;
            delete[] this->p4000;
            delete[] this->p5000;
            delete[] this->p6000;
            delete[] this->p7000;
            delete[] this->p8000;
            delete[] this->p9000;
            delete[] this->p10000;
            delete[] this->p11000;
            delete[] this->p12000;
            delete[] this->p13000;
            delete[] this->p14000;
            delete[] this->p15000;
            delete[] this->p16000;
            delete[] this->p17000;
            delete[] this->p18000;
            delete[] this->p19000;
            delete[] this->p20000;
            delete[] this->p21000;
            delete[] this->p22000;
            delete[] this->p23000;
            delete[] this->p24000;
            delete[] this->p25000;
            delete[] this->p26000;
            delete[] this->p27000;
            delete[] this->p28000;
            delete[] this->p29000;
            delete[] this->p30000;
            delete[] this->p31000;
            delete[] this->p32000;
            delete[] this->p33000;
            delete[] this->p34000;
            delete[] this->p35000;
            delete[] this->p36000;
            delete[] this->p37000;
            delete[] this->p38000;
            delete[] this->p39000;
            delete[] this->p40000;
            delete[] this->p41000;
            delete[] this->p42000;
            delete[] this->p43000;
            delete[] this->p44000;
            delete[] this->p45000;
            delete[] this->p46000;
            delete[] this->p47000;
            delete[] this->p48000;
            delete[] this->p49000;
            delete[] this->p50000;
            delete[] this->p51000;
            delete[] this->p52000;
            delete[] this->p53000;
            delete[] this->p54000;
            delete[] this->p55000;
            delete[] this->p56000;
            delete[] this->p57000;
            delete[] this->p58000;
            delete[] this->p59000;
            delete[] this->p60000;
            delete[] this->p61000;
            delete[] this->p62000;
            delete[] this->p63000;
            delete[] this->p64000;
            delete[] this->p65000;
            delete[] this->p66000;
        }
        void* malloc(size_t sz){
            int ts = __cilkrts_get_worker_number();
            switch((int)((sz-1) / 100)){
                case 0:  return p100 [ts].malloc();
                case 1:  return p200 [ts].malloc();
                case 2:  return p300 [ts].malloc();
                case 3:  return p400 [ts].malloc();
                case 4:  return p500 [ts].malloc();
                case 5:  return p600 [ts].malloc();
                case 6:  return p700 [ts].malloc();
                case 7:  return p800 [ts].malloc();
                case 8:  return p900 [ts].malloc();
                case 9:  return p1000[ts].malloc();
                case 10: return p1100[ts].malloc();
                case 11: return p1200[ts].malloc();
                case 12: return p1300[ts].malloc();
                case 13: return p1400[ts].malloc();
                case 14: return p1500[ts].malloc();
                case 15: return p1600[ts].malloc();
                case 16: return p1700[ts].malloc();
                case 17: return p1800[ts].malloc();
                case 18: return p1900[ts].malloc();
                case 19: return p2000[ts].malloc();
                case 27: return p2800[ts].malloc();
                case 31: return p3200[ts].malloc();
                case 63: return p6400[ts].malloc();
                default: switch((int)((sz-1) / 1000)){
                    case 2:  return p3000 [ts].malloc();
                    case 3:  return p4000 [ts].malloc();
                    case 4:  return p5000 [ts].malloc();
                    case 5:  return p6000 [ts].malloc();
                    case 6:  return p7000 [ts].malloc();
                    case 7:  return p8000 [ts].malloc();
                    case 8:  return p9000 [ts].malloc();
                    case 9:  return p10000[ts].malloc();
                    case 10: return p11000[ts].malloc();
                    case 11: return p12000[ts].malloc();
                    case 12: return p13000[ts].malloc();
                    case 13: return p14000[ts].malloc();
                    case 14: return p15000[ts].malloc();
                    case 15: return p16000[ts].malloc();
                    case 16: return p17000[ts].malloc();
                    case 17: return p18000[ts].malloc();
                    case 18: return p19000[ts].malloc();
                    case 19: return p20000[ts].malloc();
                    case 20: return p21000[ts].malloc();
                    case 21: return p22000[ts].malloc();
                    case 22: return p23000[ts].malloc();
                    case 23: return p24000[ts].malloc();
                    case 24: return p25000[ts].malloc();
                    case 25: return p26000[ts].malloc();
                    case 26: return p27000[ts].malloc();
                    case 27: return p28000[ts].malloc();
                    case 28: return p29000[ts].malloc();
                    case 29: return p30000[ts].malloc();
                    case 30: return p31000[ts].malloc();
                    case 31: return p32000[ts].malloc();
                    case 32: return p33000[ts].malloc();
                    case 33: return p34000[ts].malloc();
                    case 34: return p35000[ts].malloc();
                    case 35: return p36000[ts].malloc();
                    case 36: return p37000[ts].malloc();
                    case 37: return p38000[ts].malloc();
                    case 38: return p39000[ts].malloc();
                    case 39: return p40000[ts].malloc();
                    case 40: return p41000[ts].malloc();
                    case 41: return p42000[ts].malloc();
                    case 42: return p43000[ts].malloc();
                    case 43: return p44000[ts].malloc();
                    case 44: return p45000[ts].malloc();
                    case 45: return p46000[ts].malloc();
                    case 46: return p47000[ts].malloc();
                    case 47: return p48000[ts].malloc();
                    case 48: return p49000[ts].malloc();
                    case 49: return p50000[ts].malloc();
                    case 50: return p51000[ts].malloc();
                    case 51: return p52000[ts].malloc();
                    case 52: return p53000[ts].malloc();
                    case 53: return p54000[ts].malloc();
                    case 54: return p55000[ts].malloc();
                    case 55: return p56000[ts].malloc();
                    case 56: return p57000[ts].malloc();
                    case 57: return p58000[ts].malloc();
                    case 58: return p59000[ts].malloc();
                    case 59: return p60000[ts].malloc();
                    case 60: return p61000[ts].malloc();
                    case 61: return p62000[ts].malloc();
                    case 62: return p63000[ts].malloc();
                    case 63: return p64000[ts].malloc();
                    case 64: return p65000[ts].malloc();
                    case 65: return p66000[ts].malloc();
                    default: return std::malloc(sz);
                }
            }
        }
        void free(void* ptr, size_t sz){
            if(ptr == NULL) return;
            int ts = (int)(drand48()*(double)this->arity);
            switch((int)((sz-1) / 100)){
                case 0:  return p100 [ts].free(ptr);
                case 1:  return p200 [ts].free(ptr);
                case 2:  return p300 [ts].free(ptr);
                case 3:  return p400 [ts].free(ptr);
                case 4:  return p500 [ts].free(ptr);
                case 5:  return p600 [ts].free(ptr);
                case 6:  return p700 [ts].free(ptr);
                case 7:  return p800 [ts].free(ptr);
                case 8:  return p900 [ts].free(ptr);
                case 9:  return p1000[ts].free(ptr);
                case 10: return p1100[ts].free(ptr);
                case 11: return p1200[ts].free(ptr);
                case 12: return p1300[ts].free(ptr);
                case 13: return p1400[ts].free(ptr);
                case 14: return p1500[ts].free(ptr);
                case 15: return p1600[ts].free(ptr);
                case 16: return p1700[ts].free(ptr);
                case 17: return p1800[ts].free(ptr);
                case 18: return p1900[ts].free(ptr);
                case 19: return p2000[ts].free(ptr);
                case 27: return p2800[ts].free(ptr);
                case 31: return p3200[ts].free(ptr);
                case 63: return p6400[ts].free(ptr);
                default: switch((int)((sz-1) / 1000)){
                    case 2:  return p3000 [ts].free(ptr);
                    case 3:  return p4000 [ts].free(ptr);
                    case 4:  return p5000 [ts].free(ptr);
                    case 5:  return p6000 [ts].free(ptr);
                    case 6:  return p7000 [ts].free(ptr);
                    case 7:  return p8000 [ts].free(ptr);
                    case 8:  return p9000 [ts].free(ptr);
                    case 9:  return p10000[ts].free(ptr);
                    case 10: return p11000[ts].free(ptr);
                    case 11: return p12000[ts].free(ptr);
                    case 12: return p13000[ts].free(ptr);
                    case 13: return p14000[ts].free(ptr);
                    case 14: return p15000[ts].free(ptr);
                    case 15: return p16000[ts].free(ptr);
                    case 16: return p17000[ts].free(ptr);
                    case 17: return p18000[ts].free(ptr);
                    case 18: return p19000[ts].free(ptr);
                    case 19: return p20000[ts].free(ptr);
                    case 20: return p21000[ts].free(ptr);
                    case 21: return p22000[ts].free(ptr);
                    case 22: return p23000[ts].free(ptr);
                    case 23: return p24000[ts].free(ptr);
                    case 24: return p25000[ts].free(ptr);
                    case 25: return p26000[ts].free(ptr);
                    case 26: return p27000[ts].free(ptr);
                    case 27: return p28000[ts].free(ptr);
                    case 28: return p29000[ts].free(ptr);
                    case 29: return p30000[ts].free(ptr);
                    case 30: return p31000[ts].free(ptr);
                    case 31: return p32000[ts].free(ptr);
                    case 32: return p33000[ts].free(ptr);
                    case 33: return p34000[ts].free(ptr);
                    case 34: return p35000[ts].free(ptr);
                    case 35: return p36000[ts].free(ptr);
                    case 36: return p37000[ts].free(ptr);
                    case 37: return p38000[ts].free(ptr);
                    case 38: return p39000[ts].free(ptr);
                    case 39: return p40000[ts].free(ptr);
                    case 40: return p41000[ts].free(ptr);
                    case 41: return p42000[ts].free(ptr);
                    case 42: return p43000[ts].free(ptr);
                    case 43: return p44000[ts].free(ptr);
                    case 44: return p45000[ts].free(ptr);
                    case 45: return p46000[ts].free(ptr);
                    case 46: return p47000[ts].free(ptr);
                    case 47: return p48000[ts].free(ptr);
                    case 48: return p49000[ts].free(ptr);
                    case 49: return p50000[ts].free(ptr);
                    case 50: return p51000[ts].free(ptr);
                    case 51: return p52000[ts].free(ptr);
                    case 52: return p53000[ts].free(ptr);
                    case 53: return p54000[ts].free(ptr);
                    case 54: return p55000[ts].free(ptr);
                    case 55: return p56000[ts].free(ptr);
                    case 56: return p57000[ts].free(ptr);
                    case 57: return p58000[ts].free(ptr);
                    case 58: return p59000[ts].free(ptr);
                    case 59: return p60000[ts].free(ptr);
                    case 60: return p61000[ts].free(ptr);
                    case 61: return p62000[ts].free(ptr);
                    case 62: return p63000[ts].free(ptr);
                    case 63: return p64000[ts].free(ptr);
                    case 64: return p65000[ts].free(ptr);
                    case 65: return p66000[ts].free(ptr);
                    default: return std::free(ptr);
                }
            }
        }
    private:
        int arity;
    };

} }

namespace ambient {
    extern utils::bulk_memory& bulk_pool;
    extern utils::range_memory& range_pool;
    using utils::static_memory;
}

#endif
