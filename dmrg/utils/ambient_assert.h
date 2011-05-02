/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Alexandr Kosenkov <alexandr.kosenkov@unige.ch>
 *
 *****************************************************************************/

#ifndef AMBIENT_ASSERT_H
#define AMBIENT_ASSERT_H

void ambient_assert(int expression){
#ifdef MPI_PARALLEL
    assert(expression);
#endif
}

#endif
