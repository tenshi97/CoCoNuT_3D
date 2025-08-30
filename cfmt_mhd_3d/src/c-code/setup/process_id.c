#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

void get_process_id(int *process_id, int *pprocess_id) {


  int id;

  id = getpid();
  *process_id = id ;

  id = getppid();
  *pprocess_id = id ;

/*   printf("My pid %d \n",*process_id); */
/*   printf("My ppid %d \n",*pprocess_id); */


}
