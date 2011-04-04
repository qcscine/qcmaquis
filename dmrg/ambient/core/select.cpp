#include <ctype.h>
#include "ambient/core/select.h"
#include "ambient/core/scope_context.h"
#include "ambient/groups/group.h"
#include "ambient/core/operation/operation.h"
#include "utils/sqlite3.c"

namespace ambient {

    void scope_select(const char* sql)
    {
        groups::group* grp;
        int i, token_len, token_t;
        char* quantity;
        char* group; 
        char* as;
        char* field;
        char* master_token;
        char* breakdown_token;
        int master = 0;

        printf("SQL: %s\n", sql);

        i = sqlite3GetToken((const unsigned char*)sql, &token_t);
        if(token_t == TK_ILLEGAL) return;

        i += parseout<TK_ID>(sql, &group);
        i += parseout<TK_ID>(&sql[i], &as);
        if(as == NULL) as = (char*)"tmp";

        do{
            i += parseout<TK_ID>(&sql[i], &field);
            if(field != NULL)
            if(strcmp(field,"breakdown") == 0){ 
                i += parseout<TK_INTEGER>(&sql[i], &breakdown_token);
                int group_id = (int)strtol(breakdown_token, NULL, 10);
                if(group_id != 0){
                    i += parseout<TK_INTEGER>(&sql[i], &breakdown_token);
                    int id = (int)strtol(breakdown_token, NULL, 10);
// now just need to locate those ranks inside scope of profile
                }
            }else if(strcmp(field,"master") == 0){ 
                i += parseout<TK_INTEGER>(&sql[i], &master_token);
                master = (int)strtol(master_token, NULL, 10);
            }
        }while(field != NULL);

        grp = new groups::group(as, master, groups::group_map(group));

        if(token_t == TK_STAR){ 
            grp->add_every(1); 
        }else if(token_t == TK_FLOAT){ 
            float part = strtof(sql, NULL);
            grp->add_every((int)(1/part)); 
        }else if(token_t == TK_INTEGER){ 
            int count = (int)strtol(sql, NULL, 10);
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
