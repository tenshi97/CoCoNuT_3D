#ifdef CHECK_THREAD_AFFINITY

#define  _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <sched.h>

void get_thread_affinity(int *cpu_id) {

#ifndef IBM
  *cpu_id = sched_getcpu();
#else
  *cpu_id = 9999;
#endif

}

void get_process_affinity(int *cpu_id) {

#ifndef IBM
  cpu_set_t set;
#endif
  int ret, i;
  int cpu;

  *cpu_id = 9999 ;

#ifndef IBM
  ret = sched_getaffinity(0, sizeof(cpu_set_t), &set);

    for (i=0; i < CPU_SETSIZE; i++)
      {
      
	cpu = CPU_ISSET(i, &set);

	if (cpu == 1) { *cpu_id = i; }

      }

#endif
}

#endif
