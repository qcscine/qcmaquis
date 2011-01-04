#include <iostream>
#include "ambient/ambient.h"

using namespace std;

scheduler* scheduler::singleton = NULL;
scheduler* scheduler::instance(){
    if(!singleton) singleton = new scheduler();
    return singleton;
}
scheduler::scheduler(){};

scheduler & scheduler::operator>>(dim3 dim_distr) 
{
    this->dim_distr = dim_distr;
    this->dim_cpu = NULL;
    this->dim_gpu = NULL;
    return *this;
}
scheduler & scheduler::operator,(dim3 dim) 
{
    if(this->dim_cpu == NULL){
        this->dim_cpu = dim;
    }else if(this->dim_gpu == NULL){
        this->dim_gpu = dim;
    }
    return *this;
}

scheduler& operator>>(scheduler* instance, dim3 dim_distr) {
    return *instance >> dim_distr;
}

int main(int argc, char **argv)
{
    scheduler::instance() >> dim3(18,5), dim3(19,6), dim3(20,7);
    cout << "Settings were filled" << endl;
    return 0;
}

