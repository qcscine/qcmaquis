#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include "ambient/interface/select.h"
#include "ambient/core/smp.h"
#include "ambient/groups/group.h"
#include "utils/sqlite3.c"

namespace ambient {

    void select(const char* sql)
    {
        groups::group* grp;
        int i, token_len, token_t;
        char* group; 
        char* as;
        float part;
        int count;

        i = sqlite3GetToken((const unsigned char*)sql, &token_t);
        if(token_t == TK_ILLEGAL) return;

        i += parseout_id(sql, &group);
        i += parseout_id(&sql[i], &as);
        if(as == NULL) as = (char*)"tmp";
        grp = new groups::group(as, 0, group);

        if(token_t == TK_STAR){ 
            printf("selecting everything of %s as %s\n", group, as);
            grp->add_every(1); 
        }else if(token_t == TK_FLOAT){ 
            part = strtof(sql, NULL);
            printf("selecting %.2f of %s as %s\n", part, group, as);
            grp->add_every((int)(1/part)); 
        }else if(token_t == TK_INTEGER){ 
            count = (int)strtol(sql, NULL, 10);
            grp->add_range(0, count);
            printf("selecting asmps %d of %s as %s\n", count, group, as);
        }
        grp->commit();
        asmp.set_smp_group(as);
    }

    int parseout_id(const char* sql, char** id)
    {
        int i = 0;
        int token_len;
        int token_t;
        do{
            token_len = sqlite3GetToken((const unsigned char*)&sql[i], &token_t);
            i += token_len;
        }while(token_t != TK_ID && token_t != TK_ILLEGAL);

        if(token_t == TK_ILLEGAL) *id = NULL;
        else{
            *id = (char*)malloc(sizeof(char)*token_len);
            memcpy(*id, &sql[i-token_len], token_len*sizeof(char));
        }
        return i;
    }

}
