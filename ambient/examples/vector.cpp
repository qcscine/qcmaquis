#include "ambient/ambient.hpp"
#include "ambient/container/vector.hpp"

int main(){ using namespace ambient;

    vector<int> a(10, 13);
    vector<int> b(10, 10);
                                 { actor select(ambient::scope::begin()+1);
    a += b;
                                 }
                                 { actor select(ambient::scope::begin()); 
    a += b;
                                 }
    for(int i = 0; i < 10; i++)
    cout << "After sync: " << load(a).data[i] << ";\n";

    cout << "Messing with lambda:\n";
    int num = 13;

    ambient::for_each( a.begin(), a.end(), [&] (int& val){ val += num; } );
    ambient::for_each( a.begin(), a.end(), [&] (int& val){ val += num; } );

    ambient::transform( a.begin(), a.end(), b.begin(), a.begin(), [](double a, double b){ return a + b; } );

    ambient::async([&](const vector<int>& val, const vector<int>& val2){ 
            size_t size = get_length(val)-1;
            int* val_ = versioned(val).data;
            int* val2_ = versioned(val2).data;
            for(size_t i = 0; i < size; ++i) std::cout << val_[i] << " vs " << val2_[i] << "; num is " << num << "\n";
    }, a, b);

    ambient::async([&](const vector<int>& val){ 
            size_t size = get_length(val)-1;
            int* val_ = versioned(val).data;
            for(size_t i = 0; i < size; ++i) std::cout << val_[i] << "\n";
    }, a);


    ambient::sync();


    std::vector<vector<int>* > list; list.reserve(100);
    for(int i = 0; i < 100; i++){
        actor select(ambient::scope::begin()+i % 2);
        list.push_back(new vector<int>(100+i, 13));
    }

    int delta = 11;

    ambient::threaded_for_each(0, (int)list.size(), [&](int i){
        actor select(ambient::scope::balance(i, 100));
        for(int k = 0; k < 1000; k++)
        ambient::for_each( list[i]->begin(), list[i]->end(), [&] (int& val){ val += delta; } );
    });

    for(int i = 0; i < 100; i++){
        actor select(ambient::scope::begin());
        for(int k = 0; k < 1000; k++)
        ambient::for_each( list[i]->begin(), list[i]->end(), [&] (int& val){ val += delta; } );
    }

    for(int i = 0; i < 100; i++){
        ambient::async([&](const vector<int>& val){ 
            std::cout << versioned(val).data[0] << "\n";
        }, *list[i]);
        ambient::sync();
    }

    for(int i = 0; i < 100; i++) delete list[i];

    /*// possible todo:
    ambient::for_(a.begin(), b.begin(), [&](parallel_iterator& a_, parallel_iterator& b_)
    { 
        for(int i = 0; i < get_length(a); i++){
            a_[i] += b_[i];
        }
    });*/

    return 0;
}

