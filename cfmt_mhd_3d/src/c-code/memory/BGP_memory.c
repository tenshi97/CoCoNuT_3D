#if defined(CHECK_MEMORY) && defined(BLUEGENE)
#include <sys/resource.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#include <spi/kernel_interface.h>

static _BGP_Personality_t mybgp;

unsigned  bg_coreMB() {
  unsigned procMB, coreMB;

  Kernel_GetPersonality(&mybgp, sizeof(_BGP_Personality_t));
  
  procMB = BGP_Personality_DDRSizeMB(&mybgp);

  coreMB = procMB/Kernel_ProcessCount();

  return coreMB;
  
  
}

unsigned bg_stacksizeMB() {
  
  unsigned int memory_size = 0 ;
  
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &memory_size);

  memory_size = memory_size/1024/1024;
  
  return memory_size;
 
}


unsigned bg_stackavailMB() {
  
  unsigned int memory_size = 0 ;

  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &memory_size);

  memory_size = memory_size/1024/1024;
  
  return memory_size;
  
}


#endif /* definded(CHECK_MEMORY) && defined(BLUEGENE) */
