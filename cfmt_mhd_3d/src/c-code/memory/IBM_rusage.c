//#ifdef CHECK_MEMORY
#if defined(IBM) || defined(BLUEGENE)

#include <sys/resource.h>
#include <unistd.h>
#include <stdio.h>
#include <malloc.h>

#ifdef BLUEGENE
#include <bgp_SPI.h>
extern unsigned bg_coreMB();
extern unsigned bg_stacksizeMB();
extern unsigned  bg_stackavailMB();
#endif

void get_rusage_info(int *max_mem, int *heap_mem, int *heap_avail, int *heap_max, int *stack_mem, int *stack_avail, int *mempercore) {

#ifdef BLUEGENE
  uint32_t allocated_memory = 0;
  unsigned int coreMB;
#endif


  struct rusage RU;
  unsigned int top_of_heap;
  unsigned int stack;
  unsigned int tmp;


  int nullify = 0;
  

  /* Get the maximum resident size of the process */
  getrusage(RUSAGE_SELF,&RU);
  tmp=RU.ru_maxrss;
  /* on IBM Bluegene this is is kB; we want Mb*/
  tmp=tmp/1024;
  *max_mem=tmp;


#ifndef BLUEGENE
  /* on aix this is already in MB */
  top_of_heap = top_of_heap ; /* /1024/1024 ; */
  *heap_mem   = (int) top_of_heap;
  *heap_avail = nullify;
  *heap_max   = nullify;
  
  tmp        = (unsigned int) &stack /1024/1024;
  *stack_mem = (int) tmp;
  
#else

  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &allocated_memory);
  *heap_mem        = (int) allocated_memory /1024/1024;

  allocated_memory = 0;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &allocated_memory);
  *heap_avail      = allocated_memory/1024/1024;

  allocated_memory = 0;
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPMAX, &allocated_memory);
  *heap_max        = allocated_memory/1024/1024;

 
  coreMB =  bg_coreMB(); 

  *mempercore      = (int) coreMB;
  
  stack = bg_stacksizeMB();
  *stack_mem  = (int) stack;

  stack = 0;
  stack = bg_stackavailMB();
  *stack_avail = (int) stack;

#endif /* BLUEGENE */


}
#endif /* IBM */
//#endif /* CHECK_MEMORY */
