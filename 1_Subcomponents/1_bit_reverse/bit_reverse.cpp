/*
The bit_reverse part of FFT.
INPUT:
	In[]:

OUTPUT:
	Out[]:
 */

#include "bit_reverse.h"
#include <stdio.h>

unsigned int reverse_bits(unsigned int input)
{
	int i, rev = 0;
	for (i = 0; i < M; i++) {
		#pragma HLS unroll
		rev = (rev << 1) | (input & 1);
		input = input >> 1;
	}
	return rev;
}

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE]){

	unsigned int reversed;
	int x;
	DTYPE temp;
	for(x = 0; x < SIZE; x++)
	{
		reversed = reverse_bits(x);
		#pragma HLS PIPELINE II=1
		if(x < reversed)
		{
			//#pragma HLS pipeline
			temp = X_R[x];
			X_R[x] = X_R[reversed];
			X_R[reversed] = temp;

			temp = X_I[x];
			X_I[x] = X_I[reversed];
			X_I[reversed] = temp;
		}
	}

}



