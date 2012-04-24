#ifndef AMBIENT_UTILS_CEIL
#define AMBIENT_UTILS_CEIL
#define __a_ceil(x) (((double)x-(int)x) == 0 ? (int)x : (int)x+1)
#endif
