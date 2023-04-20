#include <iostream>
#include <limits>

extern "C"
{
  double ddot_ (const int64_t*, const double*, const int64_t*, const double*, const int64_t*);
}

int main (){
  int64_t inc = int64_t(3) + int64_t(2) * std::numeric_limits<int32_t>::max();
  int64_t n = 3;
  double* a = new double[3];
  a[0] = 1.0;
  a[1] = 1.0;
  a[2] = 1.0;
  // a 32 bit interface will reinterpret inc to 1 and work,
  // a 64 bit interface will take the number as is and fail with a segfault
  double d = ddot_(&n,a,&inc,a,&inc);
  delete[] a;
  std::cout << "inc = " << inc << std::endl;
  std::cout << "d   = " << d << std::endl;
  return 1;
}
