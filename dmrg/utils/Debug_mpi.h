/*
 *  Debug_mpi.h
 *  ambient_xcode
 *
 *  Created by Tim Ewart on 14.03.11.
 *  Copyright 2011 University of Geneva. All rights reserved.
 *
 */

#ifndef __DEBUG_MPI__
#define __DEBUG_MPI__

/*
 How to use :
 1- place a break point in the code where you need information
 2- start the code with mpi, you get the pid of the process, attach gdb to the process (Xcode->Run->attach to process)
 3- set the state variable to 1 (gdb command : set var  state = 1)

 Enjoy the graphic debuger, unfortunately only on one process in xcode, the second one must be controled
 by hand using gdb/terminal
 
*/

class breakpoint
{
public:
	breakpoint():state(0)
	{
		gethostname(hostname, sizeof(hostname));
		printf("PID %d on %s ready for attach\n", getpid(), hostname);
		fflush(stdout);
		stop();
	}
	
	void stop()
	{
		while(state == 0)
			sleep(1);
	}

private:	
	int state;
	char hostname[256];
	
};


#endif 

