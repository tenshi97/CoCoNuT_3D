#include <string.h>
#include <stdio.h>
#include <stdlib.h>

int test_byte_order(void)
{
#define ENDIANBIG      0
#define ENDIANLITTLE   1

   short int word = 0x0001;
   char *byte = (char *) &word;
   return(byte[0] ? ENDIANLITTLE : ENDIANBIG);

   /* return BIG if system works with big endian, otherwise LITTLE */

}
/*-----------------------------------------------------------------------------*/
/*---------------------- Routine to swap ENDIAN -------------------------------*/
/*-------- char *data:    Pointer to the data ---------------------------------*/
/*-------- int n:         Number of elements to swap --------------------------*/
/*-------- int m:         Size of single element to swap ----------------------*/
/*--------                int,float = 4 ---------------------------------------*/
/*--------                double    = 8 ---------------------------------------*/
/*-----------------------------------------------------------------------------*/
void swap_nbyte(char *data, int n, int m)
{


   int i, j;
   char old_data[16];
   
   for(j = 0; j < n; j++)
     {
       memcpy(&old_data[0], &data[j * m], m);
       for(i = 0; i < m; i++)
	 {
	   data[j * m + i] = old_data[m - i - 1];
	 }
     }
     
}
