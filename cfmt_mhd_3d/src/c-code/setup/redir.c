#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

#define NAME_LENGTH 4096
#define FILENAME "./mpi_stdout/std%3s_rank%05d.txt"

FILE *tout, *terr;
void dup_filename(char *filename, int dupfd);
void dup_fd(int fd, int dupfd);


void redirect_stdout(int *myproc) {
  char buf[NAME_LENGTH];

  if (*myproc == 0) {
#ifndef BLUEGENE
    snprintf(buf, NAME_LENGTH, "tee " FILENAME, "out", *myproc);
    tout = popen(buf, "w");
    dup_fd(fileno(tout), 1);

    snprintf(buf, NAME_LENGTH, "tee " FILENAME, "err", *myproc);
    terr = popen(buf, "w");
    dup_fd(fileno(terr), 2);
#endif /* BLUEGENE */
  } else {
    snprintf(buf, NAME_LENGTH, FILENAME, "out", *myproc);
    dup_filename(buf, 1);

    snprintf(buf, NAME_LENGTH, FILENAME, "err", *myproc);
    dup_filename(buf, 2);
  }

  return;
}

/* Redirect file descriptor dupfd to file filename */
void dup_filename(char *filename, int dupfd) {
  int fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if(fd < 0) {
    perror("open()");
    exit(1);
  }
  dup_fd(fd, dupfd);
}

/* Redirect file descriptor dupfd to file descriptor fd */
void dup_fd(int fd, int dupfd) {
  if(dup2(fd,dupfd) < 0) {
    perror("dup2()");
    exit(1);
  }
}
