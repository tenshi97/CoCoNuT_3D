//#ifdef CHECK_MEMORY
#ifndef IBM

#include <unistd.h>
#include <stdio.h>
#include <malloc.h>
#include <sys/resource.h>

void print_malloc_stat(const char *name) {
    fprintf(stderr, "Malloc_stats for %s \n",name);
    malloc_stats();

}
#endif
//#endif /* CHECK_MEMORY */
