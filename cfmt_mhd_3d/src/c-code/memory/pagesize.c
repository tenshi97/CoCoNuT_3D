//#ifdef CHECK_MEMORY
#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <malloc.h>

void get_pagesize(int *pagesize) {
  
  
  *pagesize = sysconf(_SC_PAGESIZE);
  

/* printf("Pagesize on this system %d \n\n", *pagesize ); */

}
//#endif /* CHECK_MEMORY */
