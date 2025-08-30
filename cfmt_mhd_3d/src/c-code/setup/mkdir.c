#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>

int _mkdirifnotexists(const char *dir) {
    struct stat s;
    if (stat(dir, &s) != 0) {
        if (errno == ENOENT) {
            if (mkdir(dir, 0755) != 0) {
                perror("mkdir");
                return 0;
            } else {
                return 1;
            }
        } else {
            perror("stat()");
	    return 0;
        }
    } else if (!S_ISDIR(s.st_mode)) {
        fprintf(stderr, "\"%s\" does exist and is not a directory\n", dir);
        return 0;
    } else {
        return 1;
    }
}

int create_directories(void) {
#ifdef MPI_HYDRO
    if (!_mkdirifnotexists("mpi_stdout")) return 0;
#endif
    if (!_mkdirifnotexists("output")) return 0;
    if (!_mkdirifnotexists("restart")) return 0;
    if (!_mkdirifnotexists("diagnostics")) return 0;
    return 1;
}
