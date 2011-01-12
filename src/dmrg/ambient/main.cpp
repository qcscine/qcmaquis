#include <iostream>

#include "ambient/ambient.h"

using namespace std;

int main(int argc, char **argv)
{
    scheduler::instance()->initialize();
    scheduler::instance() >> dim3(18,5), dim3(19,6), dim3(20,7);
    scheduler::instance()->finalize();
    return 0;
}
