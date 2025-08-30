#ifdef WITH_CALL_TRACE

#include <stdio.h>
#include <assert.h>
#include <execinfo.h>
#ifdef WITH_OPENMP
#include <omp.h>
#endif

const char* __attribute__ ((__no_instrument_function__)) basename(const char *string) {
  int last_slash = -1, i;
  for(i = 0; string[i] != 0; i++) {
    if (string[i] == '/') {
      last_slash = i;
    }
  }
  return &string[last_slash+1];
}

int recursion_level = 0;
# pragma omp threadprivate (recursion_level)

void __attribute__ ((__no_instrument_function__))__cyg_profile_func_enter (void *this_fn, void *call_site) {
  int thread_num = 0, num_threads = 1;
  char end[64];
  const char *func, *file; unsigned line;
  const char *call_func, *call_file; unsigned call_line;
  char **syms;
  void *stack[2];

#ifdef WITH_OPENMP
  thread_num = omp_get_thread_num();
  num_threads = omp_get_num_threads();
#endif

  recursion_level++;

  if(num_threads > 1) {
    snprintf(end, sizeof(end), " @ %d/%d\n", thread_num, num_threads);
  } else {
    snprintf(end, sizeof(end), "\n");
  }

  stack[0] = this_fn;
  stack[1] = call_site;
  syms = backtrace_symbols(stack, 2);

  fprintf(stderr, "%*scall %s from %s%s", recursion_level * 2, "", syms[0], syms[1], end);

  free(syms);
}

void __attribute__ ((__no_instrument_function__))__cyg_profile_func_exit (void *this_fn, void *call_site) {
  int thread_num = 0, num_threads = 1;
  char end[64];
  char **syms;
  void *stack[1];
#ifdef WITH_OPENMP
  thread_num = omp_get_thread_num();
  num_threads = omp_get_num_threads();
#endif

  if(num_threads > 1) {
    snprintf(end, sizeof(end), " @ %d/%d\n", thread_num, num_threads);
  } else {
    snprintf(end, sizeof(end), "\n");
  }

  stack[0] = this_fn;
  syms = backtrace_symbols(stack, 1);

  fprintf(stderr, "%*sreturn %s%s", recursion_level * 2, "", syms[0], end);

  free(syms);

  recursion_level--;
}

#endif
