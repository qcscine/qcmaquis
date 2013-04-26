//
//  main.cpp
//  TOTO
//
//  Created by Tim Ewart on 13/11/12.
//  Copyright (c) 2012 Tim Ewart. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <vector>
using std::vector;


    struct single_coefficient_task {
        unsigned int step_count; // how many time the coeff wil be calculate
        unsigned int output_degree_x;
        unsigned int output_degree_y;
        unsigned int output_degree_z;
        unsigned int output_degree_w;
    };


        int CalculateStepCount(int output_degree_x,int output_degree_y, int Order){
             int sum(0);
            for(int i1 = std::max(0,(output_degree_x-Order)); i1 <= std::min(output_degree_x, Order); ++i1)
                sum += std::min(output_degree_y,(Order-i1)) - std::max(0,(output_degree_x+output_degree_y-Order-i1))+1;
            return sum;
            }        
       
        //3 variables
        int CalculateStepCount(int output_degree_x,int output_degree_y, int output_degree_z, int Order){
                 int sum(0);
                for(int i1 = std::max(0,(output_degree_x-Order)); i1 <= std::min(output_degree_x, Order); ++i1)
                    for(int i2 = std::max(0,(output_degree_x+output_degree_y-Order)-i1); i2 <= std::min(output_degree_y, Order-i1); ++i2)
                        sum += std::min(output_degree_z,(Order-i1-i2)) - std::max(0,(output_degree_x+output_degree_y+output_degree_z-Order-i1-i2))+1;
                return sum;
            }        

        //4 variables
        int CalculateStepCount(int output_degree_x,int output_degree_y, int output_degree_z, int output_degree_w, int Order){
                 int sum(0);
                for(int i1 = std::max(0,(output_degree_x-Order)); i1 <= std::min(output_degree_x, Order); ++i1)
                    for(int i2 = std::max(0,(output_degree_x+output_degree_y-Order)-i1); i2 <= std::min(output_degree_y, Order-i1); ++i2)
                        for(int i3 = std::max(0,(output_degree_x+output_degree_y+output_degree_z-Order)-i1-i2); i3 <= std::min(output_degree_z, Order-i1-i2); ++i3)
                            sum += std::min(output_degree_w,(Order-i1-i2-i3)) - std::max(0,(output_degree_x+output_degree_y+output_degree_z+output_degree_w-Order-i1-i2-i3))+1;
                return sum;
            }    


        template<int NumVars>
        struct BuildTaskList_helper;

        template<>
        struct BuildTaskList_helper<1>{
            static void BuildTaskList(std::vector<single_coefficient_task > & VecCoeff, int Order){
                for( int degree_x = 0; degree_x <(2*Order+1); ++degree_x){
                    single_coefficient_task& task = VecCoeff[degree_x];
                    task.output_degree_x = degree_x;
                    task.output_degree_y = 0;
                    task.output_degree_z = 0;
                    task.output_degree_w = 0;
                    task.step_count =   (std::min< int>(((2*Order+1) - 1) - degree_x, degree_x) + 1);
                }
             }
        };
        // two variables 'x','y'
        template<>
        struct BuildTaskList_helper<2>{
            static void BuildTaskList(std::vector<single_coefficient_task > & VecCoeff, int Order){
                for( int degree_x = 0; degree_x <(2*Order+1); ++degree_x)
                    for( int degree_y = 0; degree_y <(2*Order+1) - degree_x; ++degree_y){
                        single_coefficient_task& task = VecCoeff[(2*Order+1)*degree_x - (degree_x*degree_x-degree_x)/2 + degree_y];
                        task.output_degree_x = degree_x;
                        task.output_degree_y = degree_y;
                        task.output_degree_z = 0;
                        task.output_degree_w = 0;
                        task.step_count = CalculateStepCount(degree_x,degree_y, Order);
                    }
            }
        };

      // tree variables 'x','y','z'
              template<>
        struct BuildTaskList_helper<3>{
            static void BuildTaskList(vector<single_coefficient_task > & VecCoeff, int Order){
                for( int degree_x = 0; degree_x <(2*Order+1); ++degree_x)
                    for( int degree_y = 0; degree_y <(2*Order+1) - degree_x; ++degree_y)
                        for( int degree_z = 0; degree_z <(2*Order+1) - degree_x - degree_y; ++degree_z){

                            single_coefficient_task& task = VecCoeff[(degree_x*(degree_x*degree_x - 3*degree_x*((2*Order+1)+1)
                                                                                  + 3*(2*Order+1)*((2*Order+1)+2) +2))/6
                                                                                  + ((2*Order+1) - degree_x)*degree_y - (degree_y*degree_y-degree_y)/2 + degree_z];
                            task.output_degree_x = degree_x;
                            task.output_degree_y = degree_y;
                            task.output_degree_z = degree_z;
                            task.output_degree_w = 0;
                            task.step_count = CalculateStepCount(degree_x,degree_y,degree_z,Order);
                    }
            }
        };
        // four variables 'x','y','z','w'
        template<>
        struct BuildTaskList_helper<4>{
            static void BuildTaskList(vector<single_coefficient_task > & VecCoeff, int Order){
                for( int degree_x = 0; degree_x <(2*Order+1); ++degree_x)
                    for( int degree_y = 0; degree_y <(2*Order+1) - degree_x; ++degree_y)
                        for( int degree_z = 0; degree_z <(2*Order+1) - degree_x - degree_y; ++degree_z)
                            for( int degree_w = 0; degree_w <(2*Order+1) - degree_x - degree_y - degree_z; ++degree_w){

                                single_coefficient_task& task = VecCoeff[(degree_x*(2*(2*Order+1)+3-degree_x)*(2*(2*Order+1)*(2*Order+1)+6*(2*Order+1)+2 +degree_x*degree_x -2*(2*Order+1)*degree_x - 3*degree_x))/24
                                                                                     + (degree_y*(degree_y*degree_y - 3*degree_y*((2*Order+1)+1-degree_x) + 3*((2*Order+1)-degree_x)*((2*Order+1)+2-degree_x)+2))/6
                                                                                     +((2*Order+1)-degree_x-degree_y)*degree_z - (degree_z*degree_z-degree_z)/2 + degree_w];
                                task.output_degree_x = degree_x;
                                task.output_degree_y = degree_y;
                                task.output_degree_z = degree_z;
                                task.output_degree_w = degree_w;
                                task.step_count = CalculateStepCount(degree_x,degree_y,degree_z,degree_w, Order);
                            }
            }
        };


   template <unsigned int N>
    struct factorial {
        static unsigned int const value = N*factorial<N-1>::value;
    };

    template <>
    struct factorial<0> {
        static unsigned int const value = 1;
    };


    template<unsigned int NK, unsigned int K>
    struct size_helper {
        static unsigned int const value = NK*size_helper<NK-1,K>::value;
    };

    template <unsigned int K>
    struct size_helper<K,K> {
        static unsigned int const value = 1;
    };

    template <unsigned int N, unsigned int K>
    struct size {
        // N variables, std::max order K -> n+k-1 over k  = (n+k-1)! / ( (n-1)! k! ) combinations
        // Assustd::ming N > 0
        static unsigned int const value = size_helper<N+K-1,K>::value/factorial<N-1>::value;
    };


#define ONE 1
#define TWO 2
#define THREE 3
#define FOUR 4
double fact( int n )
{
    return n > 1 ? n * fact( n - 1 ) : 1.0 ;
}

template< int VAR>
int SIZE(int Order){
    return fact(2*Order+VAR+1-1)/(fact(VAR)*(fact(2*Order)));
}


template<int VAR>
long foo(int Order, int sizevec){
  long qqqqq = SIZE<VAR>(Order);
   vector<single_coefficient_task > tasks(SIZE<VAR>(Order));
   BuildTaskList_helper<VAR>::BuildTaskList(tasks,Order);
   std::vector<single_coefficient_task>::iterator it = tasks.begin();
   long sum(0);
   for(; it < tasks.end(); ++it)
       sum += (*it).step_count;
    sum *= 32;
   sum *= sizevec;
//   sum += (sizevec-1)*qqqqq;
    
   return sum;
}


#define ORDER 10

int main(int argc, const char * argv[])
{
/*
   int toto = size<THREE+1,2*ORDER>::value;
   int tata =  fact(2*ORDER+THREE+1-1)/(fact(3)*(fact(2*ORDER)));

   vector<int> one,two;
   vector<single_coefficient_task > tasks(size<THREE+1,2*ORDER>::value);
   BuildTaskList_helper<THREE>::BuildTaskList(tasks,10);
   
   std::vector<single_coefficient_task>::iterator it = tasks.begin();
   int sum(0);
   for(; it < tasks.end(); ++it)
       sum += (*it).step_count;
*/
  int sizevec = 4096;
  for(int i(1); i < 15;++i)
      std::cout << foo<1>(i,sizevec) << " " << foo<2>(i,sizevec) << " " << foo<3>(i,sizevec) << " " << foo<4>(i,sizevec) << std::endl;


   return 0;
}

