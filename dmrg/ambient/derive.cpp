#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <memory.h>

class type1
{
public:
   static int data;
   int camp(){ data = 13; printf("Type1 camp!\n"); }
   template<class T>
   static int hey(){
                     printf("data is %d\n", T::data); 
                     if(T::data == 0){ T::data = 1; } 
                   }
   int parentn;

};

int type1::data = 0;

class type2: public type1
{
public:
   static int data;
   int info;
   int info2;
   type2(){
//       this->parentn = 7;
       this->info = 13;
       this->info2 = 133;
   }
   int test(){
       return 0;
   }
//   int camp(){ 
//       type1::camp();
//       this->data++;
//       printf("Type2 camp! %d\n", data); 
//   }

};

int type2::data = 0;



void test(void *memory, ...)
{
    va_list fields;
    void* buffer = malloc(200);
    void* s;

    va_start(fields, memory);
    s = va_arg(fields, void*);
//    memcpy(buffer, s, 4);

//    printf("The input is: %d\n", *((int*)buffer));
    va_end(fields);
}






int main(){

type1::hey<type1>();
type1::hey<type2>();

type2* ex = new type2();
int* nums = ((int*)ex);
int arr[] = {1, 1, 1, 165};


printf("The value of info is %d %d %d %d\n", nums[0], nums[1], nums[2], (int)sizeof(arr));

test(NULL, 123);


return 0;
}

