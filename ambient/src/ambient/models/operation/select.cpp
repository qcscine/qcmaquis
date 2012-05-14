#include <ctype.h>
#include "ambient/models/operation/select.h"
#include "ambient/controllers/context.h"
#include "ambient/channels/groups/group.h"
#include "ambient/utils/sqlite3.c"
#include "ambient/utils/timings.h"

namespace ambient {

    void ctxt_select(const char* sql)
    {
        /*
        channels::group* parent_grp;
        channels::group* grp;
        channels::group* vellum = NULL;
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
        parent_grp = channels::group_map(group);
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


        grp = new channels::group(as, master, parent_grp);

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
        catch(channels::group* original){ 
            delete grp;
            free(as);
            grp = original;
        }*/

        ctxt.set_group(channel.world());
    }

    void ctxt_retain(const char* sql)
    {
    }

    template<int TT>
    int parseout(const char* sql, char** id)
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
