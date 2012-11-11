#include <iostream>
#include <vector>
#include <algorithm>

#define ORDER 10

int main(){

std::vector<int> array((2*ORDER+1)*(2*ORDER+1));

for ( int i(0) ; i < 2*ORDER+1; ++i)
    for ( int j(0) ; j < 2*ORDER+1; ++j){
       array[j+(2*ORDER+1)*i] = (std::min((2*ORDER+1-1-i),i)+1)*(std::min((2*ORDER+1-1-j),j)+1); 
   }

std::vector<int>::iterator it = array.begin();

//for( int i=0 ; i < array.size() ; ++i)
//  std::cout << i << " " << array[i] << std::endl;


std::sort(array.begin(), array.end());
std::reverse(array.begin(), array.end());

for( int i=0 ; i < array.size() ; ++i)
  std::cout << i << " " << array[i] << std::endl;

}
