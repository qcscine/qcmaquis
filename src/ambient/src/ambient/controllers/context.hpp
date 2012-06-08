#include "ambient/channels/mpi/groups/group.h"
#include "ambient/utils/timings.hpp"

extern pthread_key_t pthread_tid;

namespace ambient { namespace controllers {


    inline context::context()
    :grp(NULL)
    { 
        this->thread_block_id = (dim2*)calloc(ambient::controller.get_num_threads(), sizeof(dim2));
    }

    inline void context::set_group(group* grp){
        this->functor->set_group(grp);
        this->grp = grp;
    }

    inline void context::set_tid(size_t value){
        void* tid = pthread_getspecific(pthread_tid);
        if(tid == NULL){
            tid = malloc(sizeof(size_t));
            pthread_setspecific(pthread_tid, tid);
        }
        *(size_t*)tid = value;
    }

    template<typename T>
    inline size_t context::get_revision_base(const iteratable<T>* o){
        return o->get_thread_revision_base();
    }

} }

namespace ambient { 

    inline void ctxt_select(const char* sql)
    {
        /*
        group* parent_grp;
        group* grp;
        group* vellum = NULL;
        int i, token_len, token_t;
        char* group           = NULL;
        char* as              = NULL;
        char* field           = NULL;
        char* master_token    = NULL;
        char* breakdown_token = NULL;
        int master = 0;

        i = sqlite3GetToken((const unsigned char*)sql, &token_t);
        if(token_t == TK_ILLEGAL) return;

        i += parseout<TK_ID>(sql, &group);
        i += parseout<TK_ID>(&sql[i], &as);
        parent_grp = group_map(group);
        free(group);

        do{
            i += parseout<TK_ID>(&sql[i], &field);
            if(field != NULL)
            if(strcmp(field,"breakdown") == 0){ 
                i += parseout<TK_INTEGER>(&sql[i], &breakdown_token);
                size_t gid = (size_t)strtol(breakdown_token, NULL, 10);
                if(gid != 0){
                    i += parseout<TK_INTEGER>(&sql[i], &breakdown_token);
                    size_t id = (size_t)strtol(breakdown_token, NULL, 10);
                    vellum = ambient::model.get_revision(&gid, 1, id)->get_placement();
                }
                free(breakdown_token);
            }else if(strcmp(field,"master") == 0){ 
                i += parseout<TK_INTEGER>(&sql[i], &master_token);
                master = (int)strtol(master_token, NULL, 10);
                free(master_token);
            }
        }while(field != NULL);
        free(field);

        vellum = ctxt.get_op()->get_vellum().get_placement(); // ignoring user preferences
        if(vellum != NULL){
           if(vellum->occupied()) vellum = NULL;
           else vellum->occupy();
        }


        grp = new group(as, master, parent_grp);

        if(token_t == TK_STAR)
        {
            if(vellum != NULL) 
                grp->add_intersection(vellum);
            else 
                grp->add_every(1); 
        }else if(token_t == TK_FLOAT)
        {
            int nth = (int)(1/strtof(sql, NULL));
            if(vellum != NULL) 
                grp->add_every_intersection(vellum, nth);
            else 
                grp->add_every(nth); 
        }else if(token_t == TK_INTEGER)
        { 
            int count = (int)strtol(sql, NULL, 10);
            if(vellum != NULL){
                grp->add_intersection(vellum, &count);
                grp->add_substraction(vellum, &count);
            }else 
                grp->add(count);
        }

        try{ grp->commit(); }
        catch(group* original){ 
            delete grp;
            free(as);
            grp = original;
        }*/

        ctxt.set_group(channel.world());
    }

    template<int TT>
    inline int parseout(const char* sql, char** id)
    {
        int i = 0;
        int token_len;
        int token_t;
        do{
            token_len = sqlite3GetToken((const unsigned char*)&sql[i], &token_t);
            i += token_len;
        }while(token_t != TT && token_t != TK_ILLEGAL);

        if(token_t == TK_ILLEGAL) *id = NULL;
        else{
            *id = (char*)realloc(*id, sizeof(char)*(token_len+1));
            memcpy(*id, &sql[i-token_len], token_len*sizeof(char));
            (*id)[token_len] = 0; // end of the string
        }
        return i;
    }

}
