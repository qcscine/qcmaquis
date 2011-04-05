#include <ctype.h>
#include "ambient/core/select.h"
#include "ambient/core/scope_context.h"
#include "ambient/groups/group.h"
#include "ambient/core/operation/operation.h"
#include "utils/sqlite3.c"

namespace ambient {

    void scope_select(const char* sql)
    {
        groups::group* parent_grp;
        groups::group* grp;
        groups::group* vellum = NULL;
        int i, token_len, token_t;
        char* quantity;
        char* group; 
        char* as;
        char* field;
        char* master_token;
        char* breakdown_token;
        int master = 0;

        i = sqlite3GetToken((const unsigned char*)sql, &token_t);
        if(token_t == TK_ILLEGAL) return;

        i += parseout<TK_ID>(sql, &group);
        i += parseout<TK_ID>(&sql[i], &as);
        parent_grp = groups::group_map(group);
        if(as == NULL) as = (char*)"tmp";

        do{
            i += parseout<TK_ID>(&sql[i], &field);
            if(field != NULL)
            if(strcmp(field,"breakdown") == 0){ 
                i += parseout<TK_INTEGER>(&sql[i], &breakdown_token);
                unsigned int group_id = (unsigned int)strtol(breakdown_token, NULL, 10);
                if(group_id != 0){
                    i += parseout<TK_INTEGER>(&sql[i], &breakdown_token);
                    unsigned int id = (unsigned int)strtol(breakdown_token, NULL, 10);
                    vellum = p_profile_map.find(&group_id, 1, id)->profile->get_scope();
                }
            }else if(strcmp(field,"master") == 0){ 
                i += parseout<TK_INTEGER>(&sql[i], &master_token);
                master = (int)strtol(master_token, NULL, 10);
            }
        }while(field != NULL);

        grp = new groups::group(as, master, parent_grp);

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
        grp->commit();
        scope.set_group(grp);
        scope.get_op()->preprocess();
    }

    void scope_retain(const char* sql)
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
            *id = (char*)malloc(sizeof(char)*(token_len+1));
            memcpy(*id, &sql[i-token_len], token_len*sizeof(char));
            (*id)[token_len] = 0; // end of the string
        }
        return i;
    }

}
