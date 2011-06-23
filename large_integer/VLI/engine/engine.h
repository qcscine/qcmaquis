//
//  Header.h
//  vli
//
//  Created by Tim Ewart on 15/06/11.
//  Copyright 2011 University of Geneva. All rights reserved.
//
#include "engine/membank.h"

typedef bank::membank< vli::vli_cpu<int> > genericbank;
typedef bank::membankcpu< vli::vli_cpu<int> > cpubank;

extern genericbank* CPUBANK; 