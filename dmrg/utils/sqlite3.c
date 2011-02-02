#define sqlite3StrNICmp     sqlite3_strnicmp
#define UpperToLower        sqlite3UpperToLower
#define charMap(X)          sqlite3UpperToLower[(unsigned char)X]
#define IdChar(C)           ((sqlite3CtypeMap[(unsigned char)C]&0x46)!=0)
#define sqlite3Isspace(x)   isspace((unsigned char)(x))
#define sqlite3Isdigit(x)   isdigit((unsigned char)(x))

#define TK_SEMI                            1
#define TK_EXPLAIN                         2
#define TK_QUERY                           3
#define TK_PLAN                            4
#define TK_BEGIN                           5
#define TK_TRANSACTION                     6
#define TK_DEFERRED                        7
#define TK_IMMEDIATE                       8
#define TK_EXCLUSIVE                       9
#define TK_COMMIT                         10
#define TK_END                            11
#define TK_ROLLBACK                       12
#define TK_SAVEPOINT                      13
#define TK_RELEASE                        14
#define TK_TO                             15
#define TK_TABLE                          16
#define TK_CREATE                         17
#define TK_IF                             18
#define TK_NOT                            19
#define TK_EXISTS                         20
#define TK_TEMP                           21
#define TK_LP                             22
#define TK_RP                             23
#define TK_AS                             24
#define TK_COMMA                          25
#define TK_ID                             26
#define TK_INDEXED                        27
#define TK_ABORT                          28
#define TK_ACTION                         29
#define TK_AFTER                          30
#define TK_ANALYZE                        31
#define TK_ASC                            32
#define TK_ATTACH                         33
#define TK_BEFORE                         34
#define TK_BY                             35
#define TK_CASCADE                        36
#define TK_CAST                           37
#define TK_COLUMNKW                       38
#define TK_CONFLICT                       39
#define TK_DATABASE                       40
#define TK_DESC                           41
#define TK_DETACH                         42
#define TK_EACH                           43
#define TK_FAIL                           44
#define TK_FOR                            45
#define TK_IGNORE                         46
#define TK_INITIALLY                      47
#define TK_INSTEAD                        48
#define TK_LIKE_KW                        49
#define TK_MATCH                          50
#define TK_NO                             51
#define TK_KEY                            52
#define TK_OF                             53
#define TK_OFFSET                         54
#define TK_PRAGMA                         55
#define TK_RAISE                          56
#define TK_REPLACE                        57
#define TK_RESTRICT                       58
#define TK_ROW                            59
#define TK_TRIGGER                        60
#define TK_VACUUM                         61
#define TK_VIEW                           62
#define TK_VIRTUAL                        63
#define TK_REINDEX                        64
#define TK_RENAME                         65
#define TK_CTIME_KW                       66
#define TK_ANY                            67
#define TK_OR                             68
#define TK_AND                            69
#define TK_IS                             70
#define TK_BETWEEN                        71
#define TK_IN                             72
#define TK_ISNULL                         73
#define TK_NOTNULL                        74
#define TK_NE                             75
#define TK_EQ                             76
#define TK_GT                             77
#define TK_LE                             78
#define TK_LT                             79
#define TK_GE                             80
#define TK_ESCAPE                         81
#define TK_BITAND                         82
#define TK_BITOR                          83
#define TK_LSHIFT                         84
#define TK_RSHIFT                         85
#define TK_PLUS                           86
#define TK_MINUS                          87
#define TK_STAR                           88
#define TK_SLASH                          89
#define TK_REM                            90
#define TK_CONCAT                         91
#define TK_COLLATE                        92
#define TK_BITNOT                         93
#define TK_STRING                         94
#define TK_JOIN_KW                        95
#define TK_CONSTRAINT                     96
#define TK_DEFAULT                        97
#define TK_NULL                           98
#define TK_PRIMARY                        99
#define TK_UNIQUE                         100
#define TK_CHECK                          101
#define TK_REFERENCES                     102
#define TK_AUTOINCR                       103
#define TK_ON                             104
#define TK_INSERT                         105
#define TK_DELETE                         106
#define TK_UPDATE                         107
#define TK_SET                            108
#define TK_DEFERRABLE                     109
#define TK_FOREIGN                        110
#define TK_DROP                           111
#define TK_UNION                          112
#define TK_ALL                            113
#define TK_EXCEPT                         114
#define TK_INTERSECT                      115
#define TK_SELECT                         116
#define TK_DISTINCT                       117
#define TK_DOT                            118
#define TK_FROM                           119
#define TK_JOIN                           120
#define TK_USING                          121
#define TK_ORDER                          122
#define TK_GROUP                          123
#define TK_HAVING                         124
#define TK_LIMIT                          125
#define TK_WHERE                          126
#define TK_INTO                           127
#define TK_VALUES                         128
#define TK_INTEGER                        129
#define TK_FLOAT                          130
#define TK_BLOB                           131
#define TK_REGISTER                       132
#define TK_VARIABLE                       133
#define TK_CASE                           134
#define TK_WHEN                           135
#define TK_THEN                           136
#define TK_ELSE                           137
#define TK_INDEX                          138
#define TK_ALTER                          139
#define TK_ADD                            140
#define TK_TO_TEXT                        141
#define TK_TO_BLOB                        142
#define TK_TO_NUMERIC                     143
#define TK_TO_INT                         144
#define TK_TO_REAL                        145
#define TK_ISNOT                          146
#define TK_END_OF_FILE                    147
#define TK_ILLEGAL                        148
#define TK_SPACE                          149
#define TK_UNCLOSED_STRING                150
#define TK_FUNCTION                       151
#define TK_COLUMN                         152
#define TK_AGG_FUNCTION                   153
#define TK_AGG_COLUMN                     154
#define TK_CONST_FUNC                     155
#define TK_UMINUS                         156
#define TK_UPLUS                          157

const unsigned char sqlite3UpperToLower[] = {
      0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
     18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
     36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53,
     54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 97, 98, 99,100,101,102,103,
    104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,
    122, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103,104,105,106,107,
    108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,
    126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
    144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,
    162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,
    180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,
    198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,
    216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,
    234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,
    252,253,254,255
};
const unsigned char sqlite3CtypeMap[256] = {
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  /* 00..07    ........ */
  0x00, 0x01, 0x01, 0x01, 0x01, 0x01, 0x00, 0x00,  /* 08..0f    ........ */
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  /* 10..17    ........ */
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  /* 18..1f    ........ */
  0x01, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00,  /* 20..27     !"#$%&' */
  0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  /* 28..2f    ()*+,-./ */
  0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c,  /* 30..37    01234567 */
  0x0c, 0x0c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,  /* 38..3f    89:;<=>? */

  0x00, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x0a, 0x02,  /* 40..47    @ABCDEFG */
  0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02,  /* 48..4f    HIJKLMNO */
  0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02, 0x02,  /* 50..57    PQRSTUVW */
  0x02, 0x02, 0x02, 0x00, 0x00, 0x00, 0x00, 0x40,  /* 58..5f    XYZ[\]^_ */
  0x00, 0x2a, 0x2a, 0x2a, 0x2a, 0x2a, 0x2a, 0x22,  /* 60..67    `abcdefg */
  0x22, 0x22, 0x22, 0x22, 0x22, 0x22, 0x22, 0x22,  /* 68..6f    hijklmno */
  0x22, 0x22, 0x22, 0x22, 0x22, 0x22, 0x22, 0x22,  /* 70..77    pqrstuvw */
  0x22, 0x22, 0x22, 0x00, 0x00, 0x00, 0x00, 0x00,  /* 78..7f    xyz{|}~. */

  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* 80..87    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* 88..8f    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* 90..97    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* 98..9f    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* a0..a7    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* a8..af    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* b0..b7    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* b8..bf    ........ */

  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* c0..c7    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* c8..cf    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* d0..d7    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* d8..df    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* e0..e7    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* e8..ef    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40,  /* f0..f7    ........ */
  0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40, 0x40   /* f8..ff    ........ */
};

int sqlite3_strnicmp(const char *zLeft, const char *zRight, int N){
  register unsigned char *a, *b;
  a = (unsigned char *)zLeft;
  b = (unsigned char *)zRight;
  while( N-- > 0 && *a!=0 && UpperToLower[*a]==UpperToLower[*b]){ a++; b++; }
  return N<0 ? 0 : UpperToLower[*a] - UpperToLower[*b];
}


/***** This file contains automatically generated code ******
**
** The code in this file has been automatically generated by
**
**   sqlite/tool/mkkeywordhash.c
**
** The code in this file implements a function that determines whether
** or not a given identifier is really an SQL keyword.  The same thing
** might be implemented more directly using a hand-written hash table.
** But by using this automatically generated code, the size of the code
** is substantially reduced.  This is important for embedded applications
** on platforms with limited memory.
*/
/* Hash score: 175 */
static int keywordCode(const char *z, int n){
  /* zText[] encodes 811 bytes of keywords in 541 bytes */
  /*   REINDEXEDESCAPEACHECKEYBEFOREIGNOREGEXPLAINSTEADDATABASELECT       */
  /*   ABLEFTHENDEFERRABLELSEXCEPTRANSACTIONATURALTERAISEXCLUSIVE         */
  /*   XISTSAVEPOINTERSECTRIGGEREFERENCESCONSTRAINTOFFSETEMPORARY         */
  /*   UNIQUERYATTACHAVINGROUPDATEBEGINNERELEASEBETWEENOTNULLIKE          */
  /*   CASCADELETECASECOLLATECREATECURRENT_DATEDETACHIMMEDIATEJOIN        */
  /*   SERTMATCHPLANALYZEPRAGMABORTVALUESVIRTUALIMITWHENWHERENAME         */
  /*   AFTEREPLACEANDEFAULTAUTOINCREMENTCASTCOLUMNCOMMITCONFLICTCROSS     */
  /*   CURRENT_TIMESTAMPRIMARYDEFERREDISTINCTDROPFAILFROMFULLGLOBYIF      */
  /*   ISNULLORDERESTRICTOUTERIGHTROLLBACKROWUNIONUSINGVACUUMVIEW         */
  /*   INITIALLY                                                          */
  static const char zText[540] = {
    'R','E','I','N','D','E','X','E','D','E','S','C','A','P','E','A','C','H',
    'E','C','K','E','Y','B','E','F','O','R','E','I','G','N','O','R','E','G',
    'E','X','P','L','A','I','N','S','T','E','A','D','D','A','T','A','B','A',
    'S','E','L','E','C','T','A','B','L','E','F','T','H','E','N','D','E','F',
    'E','R','R','A','B','L','E','L','S','E','X','C','E','P','T','R','A','N',
    'S','A','C','T','I','O','N','A','T','U','R','A','L','T','E','R','A','I',
    'S','E','X','C','L','U','S','I','V','E','X','I','S','T','S','A','V','E',
    'P','O','I','N','T','E','R','S','E','C','T','R','I','G','G','E','R','E',
    'F','E','R','E','N','C','E','S','C','O','N','S','T','R','A','I','N','T',
    'O','F','F','S','E','T','E','M','P','O','R','A','R','Y','U','N','I','Q',
    'U','E','R','Y','A','T','T','A','C','H','A','V','I','N','G','R','O','U',
    'P','D','A','T','E','B','E','G','I','N','N','E','R','E','L','E','A','S',
    'E','B','E','T','W','E','E','N','O','T','N','U','L','L','I','K','E','C',
    'A','S','C','A','D','E','L','E','T','E','C','A','S','E','C','O','L','L',
    'A','T','E','C','R','E','A','T','E','C','U','R','R','E','N','T','_','D',
    'A','T','E','D','E','T','A','C','H','I','M','M','E','D','I','A','T','E',
    'J','O','I','N','S','E','R','T','M','A','T','C','H','P','L','A','N','A',
    'L','Y','Z','E','P','R','A','G','M','A','B','O','R','T','V','A','L','U',
    'E','S','V','I','R','T','U','A','L','I','M','I','T','W','H','E','N','W',
    'H','E','R','E','N','A','M','E','A','F','T','E','R','E','P','L','A','C',
    'E','A','N','D','E','F','A','U','L','T','A','U','T','O','I','N','C','R',
    'E','M','E','N','T','C','A','S','T','C','O','L','U','M','N','C','O','M',
    'M','I','T','C','O','N','F','L','I','C','T','C','R','O','S','S','C','U',
    'R','R','E','N','T','_','T','I','M','E','S','T','A','M','P','R','I','M',
    'A','R','Y','D','E','F','E','R','R','E','D','I','S','T','I','N','C','T',
    'D','R','O','P','F','A','I','L','F','R','O','M','F','U','L','L','G','L',
    'O','B','Y','I','F','I','S','N','U','L','L','O','R','D','E','R','E','S',
    'T','R','I','C','T','O','U','T','E','R','I','G','H','T','R','O','L','L',
    'B','A','C','K','R','O','W','U','N','I','O','N','U','S','I','N','G','V',
    'A','C','U','U','M','V','I','E','W','I','N','I','T','I','A','L','L','Y',
  };
  static const unsigned char aHash[127] = {
      72, 101, 114,  70,   0,  45,   0,   0,  78,   0,  73,   0,   0,
      42,  12,  74,  15,   0, 113,  81,  50, 108,   0,  19,   0,   0,
     118,   0, 116, 111,   0,  22,  89,   0,   9,   0,   0,  66,  67,
       0,  65,   6,   0,  48,  86,  98,   0, 115,  97,   0,   0,  44,
       0,  99,  24,   0,  17,   0, 119,  49,  23,   0,   5, 106,  25,
      92,   0,   0, 121, 102,  56, 120,  53,  28,  51,   0,  87,   0,
      96,  26,   0,  95,   0,   0,   0,  91,  88,  93,  84, 105,  14,
      39, 104,   0,  77,   0,  18,  85, 107,  32,   0, 117,  76, 109,
      58,  46,  80,   0,   0,  90,  40,   0, 112,   0,  36,   0,   0,
      29,   0,  82,  59,  60,   0,  20,  57,   0,  52,
  };
  static const unsigned char aNext[121] = {
       0,   0,   0,   0,   4,   0,   0,   0,   0,   0,   0,   0,   0,
       0,   2,   0,   0,   0,   0,   0,   0,  13,   0,   0,   0,   0,
       0,   7,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
       0,   0,   0,   0,  33,   0,  21,   0,   0,   0,  43,   3,  47,
       0,   0,   0,   0,  30,   0,  54,   0,  38,   0,   0,   0,   1,
      62,   0,   0,  63,   0,  41,   0,   0,   0,   0,   0,   0,   0,
      61,   0,   0,   0,   0,  31,  55,  16,  34,  10,   0,   0,   0,
       0,   0,   0,   0,  11,  68,  75,   0,   8,   0, 100,  94,   0,
     103,   0,  83,   0,  71,   0,   0, 110,  27,  37,  69,  79,   0,
      35,  64,   0,   0,
  };
  static const unsigned char aLen[121] = {
       7,   7,   5,   4,   6,   4,   5,   3,   6,   7,   3,   6,   6,
       7,   7,   3,   8,   2,   6,   5,   4,   4,   3,  10,   4,   6,
      11,   6,   2,   7,   5,   5,   9,   6,   9,   9,   7,  10,  10,
       4,   6,   2,   3,   9,   4,   2,   6,   5,   6,   6,   5,   6,
       5,   5,   7,   7,   7,   3,   2,   4,   4,   7,   3,   6,   4,
       7,   6,  12,   6,   9,   4,   6,   5,   4,   7,   6,   5,   6,
       7,   5,   4,   5,   6,   5,   7,   3,   7,  13,   2,   2,   4,
       6,   6,   8,   5,  17,  12,   7,   8,   8,   2,   4,   4,   4,
       4,   4,   2,   2,   6,   5,   8,   5,   5,   8,   3,   5,   5,
       6,   4,   9,   3,
  };
  static const unsigned short int aOffset[121] = {
       0,   2,   2,   8,   9,  14,  16,  20,  23,  25,  25,  29,  33,
      36,  41,  46,  48,  53,  54,  59,  62,  65,  67,  69,  78,  81,
      86,  91,  95,  96, 101, 105, 109, 117, 122, 128, 136, 142, 152,
     159, 162, 162, 165, 167, 167, 171, 176, 179, 184, 189, 194, 197,
     203, 206, 210, 217, 223, 223, 223, 226, 229, 233, 234, 238, 244,
     248, 255, 261, 273, 279, 288, 290, 296, 301, 303, 310, 315, 320,
     326, 332, 337, 341, 344, 350, 354, 361, 363, 370, 372, 374, 383,
     387, 393, 399, 407, 412, 412, 428, 435, 442, 443, 450, 454, 458,
     462, 466, 469, 471, 473, 479, 483, 491, 495, 500, 508, 511, 516,
     521, 527, 531, 536,
  };
  static const unsigned char aCode[121] = {
    TK_REINDEX,    TK_INDEXED,    TK_INDEX,      TK_DESC,       TK_ESCAPE,     
    TK_EACH,       TK_CHECK,      TK_KEY,        TK_BEFORE,     TK_FOREIGN,    
    TK_FOR,        TK_IGNORE,     TK_LIKE_KW,    TK_EXPLAIN,    TK_INSTEAD,    
    TK_ADD,        TK_DATABASE,   TK_AS,         TK_SELECT,     TK_TABLE,      
    TK_JOIN_KW,    TK_THEN,       TK_END,        TK_DEFERRABLE, TK_ELSE,       
    TK_EXCEPT,     TK_TRANSACTION,TK_ACTION,     TK_ON,         TK_JOIN_KW,    
    TK_ALTER,      TK_RAISE,      TK_EXCLUSIVE,  TK_EXISTS,     TK_SAVEPOINT,  
    TK_INTERSECT,  TK_TRIGGER,    TK_REFERENCES, TK_CONSTRAINT, TK_INTO,       
    TK_OFFSET,     TK_OF,         TK_SET,        TK_TEMP,       TK_TEMP,       
    TK_OR,         TK_UNIQUE,     TK_QUERY,      TK_ATTACH,     TK_HAVING,     
    TK_GROUP,      TK_UPDATE,     TK_BEGIN,      TK_JOIN_KW,    TK_RELEASE,    
    TK_BETWEEN,    TK_NOTNULL,    TK_NOT,        TK_NO,         TK_NULL,       
    TK_LIKE_KW,    TK_CASCADE,    TK_ASC,        TK_DELETE,     TK_CASE,       
    TK_COLLATE,    TK_CREATE,     TK_CTIME_KW,   TK_DETACH,     TK_IMMEDIATE,  
    TK_JOIN,       TK_INSERT,     TK_MATCH,      TK_PLAN,       TK_ANALYZE,    
    TK_PRAGMA,     TK_ABORT,      TK_VALUES,     TK_VIRTUAL,    TK_LIMIT,      
    TK_WHEN,       TK_WHERE,      TK_RENAME,     TK_AFTER,      TK_REPLACE,    
    TK_AND,        TK_DEFAULT,    TK_AUTOINCR,   TK_TO,         TK_IN,         
    TK_CAST,       TK_COLUMNKW,   TK_COMMIT,     TK_CONFLICT,   TK_JOIN_KW,    
    TK_CTIME_KW,   TK_CTIME_KW,   TK_PRIMARY,    TK_DEFERRED,   TK_DISTINCT,   
    TK_IS,         TK_DROP,       TK_FAIL,       TK_FROM,       TK_JOIN_KW,    
    TK_LIKE_KW,    TK_BY,         TK_IF,         TK_ISNULL,     TK_ORDER,      
    TK_RESTRICT,   TK_JOIN_KW,    TK_JOIN_KW,    TK_ROLLBACK,   TK_ROW,        
    TK_UNION,      TK_USING,      TK_VACUUM,     TK_VIEW,       TK_INITIALLY,  
    TK_ALL,        
  };
  int h, i;
  if( n<2 ) return TK_ID;
  h = ((charMap(z[0])*4) ^
      (charMap(z[n-1])*3) ^
      n) % 127;
  for(i=((int)aHash[h])-1; i>=0; i=((int)aNext[i])-1){
    if( aLen[i]==n && sqlite3StrNICmp(&zText[aOffset[i]],z,n)==0 ){
      return aCode[i];
    }
  }
  return TK_ID;
}
/*
** Return the length of the token that begins at z[0]. 
** Store the token type in *tokenType before returning.
*/
int sqlite3GetToken(const unsigned char *z, int *tokenType){
  int i, c;
  switch( *z ){
    case ' ': case '\t': case '\n': case '\f': case '\r': {
      for(i=1; sqlite3Isspace(z[i]); i++){}
      *tokenType = TK_SPACE;
      return i;
    }
    case '-': {
      if( z[1]=='-' ){
        /* IMP: R-15891-05542 -- syntax diagram for comments */
        for(i=2; (c=z[i])!=0 && c!='\n'; i++){}
        *tokenType = TK_SPACE;   /* IMP: R-22934-25134 */
        return i;
      }
      *tokenType = TK_MINUS;
      return 1;
    }
    case '(': {
      *tokenType = TK_LP;
      return 1;
    }
    case ')': {
      *tokenType = TK_RP;
      return 1;
    }
    case ';': {
      *tokenType = TK_SEMI;
      return 1;
    }
    case '+': {
      *tokenType = TK_PLUS;
      return 1;
    }
    case '*': {
      *tokenType = TK_STAR;
      return 1;
    }
    case '/': {
      if( z[1]!='*' || z[2]==0 ){
        *tokenType = TK_SLASH;
        return 1;
      }
      /* IMP: R-15891-05542 -- syntax diagram for comments */
      for(i=3, c=z[2]; (c!='*' || z[i]!='/') && (c=z[i])!=0; i++){}
      if( c ) i++;
      *tokenType = TK_SPACE;   /* IMP: R-22934-25134 */
      return i;
    }
    case '%': {
      *tokenType = TK_REM;
      return 1;
    }
    case '=': {
      *tokenType = TK_EQ;
      return 1 + (z[1]=='=');
    }
    case '<': {
      if( (c=z[1])=='=' ){
        *tokenType = TK_LE;
        return 2;
      }else if( c=='>' ){
        *tokenType = TK_NE;
        return 2;
      }else if( c=='<' ){
        *tokenType = TK_LSHIFT;
        return 2;
      }else{
        *tokenType = TK_LT;
        return 1;
      }
    }
    case '>': {
      if( (c=z[1])=='=' ){
        *tokenType = TK_GE;
        return 2;
      }else if( c=='>' ){
        *tokenType = TK_RSHIFT;
        return 2;
      }else{
        *tokenType = TK_GT;
        return 1;
      }
    }
    case '!': {
      if( z[1]!='=' ){
        *tokenType = TK_ILLEGAL;
        return 2;
      }else{
        *tokenType = TK_NE;
        return 2;
      }
    }
    case '|': {
      if( z[1]!='|' ){
        *tokenType = TK_BITOR;
        return 1;
      }else{
        *tokenType = TK_CONCAT;
        return 2;
      }
    }
    case ',': {
      *tokenType = TK_COMMA;
      return 1;
    }
    case '&': {
      *tokenType = TK_BITAND;
      return 1;
    }
    case '~': {
      *tokenType = TK_BITNOT;
      return 1;
    }
    case '`':
    case '\'':
    case '"': {
      int delim = z[0];
      for(i=1; (c=z[i])!=0; i++){
        if( c==delim ){
          if( z[i+1]==delim ){
            i++;
          }else{
            break;
          }
        }
      }
      if( c=='\'' ){
        *tokenType = TK_STRING;
        return i+1;
      }else if( c!=0 ){
        *tokenType = TK_ID;
        return i+1;
      }else{
        *tokenType = TK_ILLEGAL;
        return i;
      }
    }
    case '.': {
      if( !sqlite3Isdigit(z[1]) )
      {
        *tokenType = TK_DOT;
        return 1;
      }
      /* If the next character is a digit, this is a floating point
      ** number that begins with ".".  Fall thru into the next case */
    }
    case '0': case '1': case '2': case '3': case '4':
    case '5': case '6': case '7': case '8': case '9': {
      *tokenType = TK_INTEGER;
      for(i=0; sqlite3Isdigit(z[i]); i++){}
      if( z[i]=='.' ){
        i++;
        while( sqlite3Isdigit(z[i]) ){ i++; }
        *tokenType = TK_FLOAT;
      }
      if( (z[i]=='e' || z[i]=='E') &&
           ( sqlite3Isdigit(z[i+1]) 
            || ((z[i+1]=='+' || z[i+1]=='-') && sqlite3Isdigit(z[i+2]))
           )
      ){
        i += 2;
        while( sqlite3Isdigit(z[i]) ){ i++; }
        *tokenType = TK_FLOAT;
      }
      while( IdChar(z[i]) ){
        *tokenType = TK_ILLEGAL;
        i++;
      }
      return i;
    }
    case '[': {
      for(i=1, c=z[0]; c!=']' && (c=z[i])!=0; i++){}
      *tokenType = c==']' ? TK_ID : TK_ILLEGAL;
      return i;
    }
    case '?': {
      *tokenType = TK_VARIABLE;
      for(i=1; sqlite3Isdigit(z[i]); i++){}
      return i;
    }
    case '#': {
      for(i=1; sqlite3Isdigit(z[i]); i++){}
      if( i>1 ){
        /* Parameters of the form #NNN (where NNN is a number) are used
        ** internally by sqlite3NestedParse.  */
        *tokenType = TK_REGISTER;
        return i;
      }
      /* Fall through into the next case if the '#' is not followed by
      ** a digit. Try to match #AAAA where AAAA is a parameter name. */
    }
    case '@':  /* For compatibility with MS SQL Server */
    case ':': {
      int n = 0;
      *tokenType = TK_VARIABLE;
      for(i=1; (c=z[i])!=0; i++){
        if( IdChar(c) ){
          n++;
        }else{
          break;
        }
      }
      if( n==0 ) *tokenType = TK_ILLEGAL;
      return i;
    }
    default: {
      if( !IdChar(*z) ){
        break;
      }
      for(i=1; IdChar(z[i]); i++){}
      *tokenType = keywordCode((char*)z, i);
      return i;
    }
  }
  *tokenType = TK_ILLEGAL;
  return 1;
}
