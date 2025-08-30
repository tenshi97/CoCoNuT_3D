#ifdef LINUX
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

void unlimit_stack(void) {
  struct rlimit rlim = {RLIM_INFINITY, RLIM_INFINITY};

  if (setrlimit(RLIMIT_STACK, &rlim) == -1) {
    perror("setrlimit()");
  }
}
#endif
