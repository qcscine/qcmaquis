#ifndef AMBIENT_CONTROLLERS_VELVET_CHAIN
#define AMBIENT_CONTROLLERS_VELVET_CHAIN

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::revision;

    class chain {
    public:
        chain(cfunctor* f){
            affinity = f->affinity;
            this->push_back(f);
        }
        void push_back(cfunctor* f){
            std::list<revision*>::iterator it; 
            this->content.push_back(f);

            std::list<revision*>& ders = f->derivatives; it = ders.begin(); 
            while(it != ders.end()) this->derivatives.push_back(*it++);

            std::list<revision*>& deps = f->dependencies; it = deps.begin(); 
            while(it != deps.end()) this->dependencies.push_back(*it++);
        }
        void execute(){
            cfunctor* op;
            std::list<cfunctor*>::iterator it = this->content.begin(); 
            while(it != this->content.end()){
                op = *it++;
                op->computation();
                std::list<revision*>& list = op->derivatives; 
                std::list<revision*>::iterator itd = list.begin(); 
                while(itd != list.end()) (*itd++)->reset_generator();
                delete op;
            }
            this->derivatives.clear();
        }
        bool ready(){
            std::list<revision*>::iterator it = this->dependencies.begin(); 
            while(it != this->dependencies.end()) if((*it++)->generator != NULL) return false;
            return true;
        }
        bool constrains(cfunctor* op){
            bool result = false;
            std::list<revision*>::iterator it = op->dependencies.begin(); 
rt:         while(it != op->dependencies.end())
            {
                std::list<revision*>::iterator itd = this->derivatives.begin(); 
                while(itd != this->derivatives.end()){
                    if(*it == *itd){
                        result = true;
                        op->dependencies.erase(it++);
                        goto rt;
                    }else
                        itd++;
                }
                it++;
            }
            return result;
        }
        int affinity;
        std::list<cfunctor*> content;
    private:
        std::list<revision*> derivatives;
        std::list<revision*> dependencies;
    };

} } }

#endif
