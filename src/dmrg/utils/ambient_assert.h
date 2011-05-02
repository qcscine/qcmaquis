/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Alexandr Kosenkov <alexandr.kosenkov@unige.ch>
 *
 *****************************************************************************/

#ifndef AMBIENT_ASSERT_H
#define AMBIENT_ASSERT_H

#ifdef MPI_PARALLEL
#define ambient_assert(...) assert( __VA_ARGS__ );
#else
#define ambient_assert(...) assert(true);
#endif

#endif
