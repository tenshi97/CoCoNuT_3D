#ifdef BACKTRACE_ON_ABORT

#include <stdio.h>
#include <execinfo.h>

void print_backtrace(void) {
    void *buffer[100];
    int n_frames;
    n_frames = backtrace(buffer, 100);
    backtrace_symbols_fd(buffer, n_frames, 1);
}

#endif /* BACKTRACE_ON_ABORT */
