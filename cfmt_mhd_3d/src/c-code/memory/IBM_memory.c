//#ifdef CHECK_MEMORY
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/resource.h>
#include <malloc.h>

void ibm_mem_free(int *mem) {
	FILE *fd;
	int pages;

	int pagesize;

	pagesize =  sysconf(_SC_PAGESIZE)/1024;

	fd = popen("vmstat -v | grep free", "r");
	if (fd == NULL) {
		perror("popen");
		*mem = 0;
		return;
	}

	fscanf(fd, "%d free pages", &pages);
	pclose(fd);


	/* convert pages to MB */
    	*mem = pages*pagesize/1024;

}
//#endif /* CHECK_MEMORY */
