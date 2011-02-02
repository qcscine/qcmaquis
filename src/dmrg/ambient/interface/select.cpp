#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include "ambient/interface/select.h"
#include "utils/sqlite3.c"

namespace ambient {

    void select(const char* sql)
    {
        int i=0;
        int token_len;
        int token_t;

        for(;;){
            token_len = sqlite3GetToken((const unsigned char*)&sql[i], &token_t);
            if(token_t == TK_ILLEGAL) break;
    // parse the result here
            printf("Length is %d and type is %d\n", token_len, token_t);
    // end of parsing
            i += token_len;
        }
    }

}
