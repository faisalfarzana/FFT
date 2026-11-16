#include"fft_stage_last.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * Function: fft_stage_first_op calculates the LAST stage of FFT.
 * 
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal

OUTPUT:
	OUT_R, OUT_I[]: Real and Imag parts of Complex signal

 */

//=======================BEGIN: FFT fft_stage_last =========================
void fft_stage_last(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

	//Write code that computes last stage of FFT
	//Insert your code here
	DTYPE temp_R;		/*temporary storage complex variable*/
	DTYPE temp_I;		/*temporary storage complex variable*/


	int i,j,k;			/* loop indexes */
	int i_lower;		/* Index of lower point in butterfly */
	int step;

	int stage = 10;
	int DFTpts;
	int numBF;			/*Butterfly Width*/

	int N2 = SIZE2;	/* N2=N>>1 */
	step=N2;
	DTYPE a, e, c, s;

		DFTpts = 1024;		// DFT = 2^stage = points in sub DFT
		numBF = DFTpts/2; 			// Butterfly WIDTHS in sub-DFT
		k=0;

		e = -6.283185307178/DFTpts;
		a = 0.0;
		// Perform butterflies for j-th stage
		butterfly:for(j=0; j<numBF; j++)
		{
			//c = cos(a);
			//s = sin(a);
			//a = a + e;

			// Compute butterflies that use same W**k
			DFTpts:for(i=j; i<SIZE; i += DFTpts)
			{
				#pragma HLS unroll
				i_lower = i + numBF;			//index of lower point in butterfly
				temp_R = X_R[i_lower]*W_real[j]- X_I[i_lower]*W_imag[j];
				temp_I = X_I[i_lower]*W_real[j]+ X_R[i_lower]*W_imag[j];

				OUT_R[i_lower] = X_R[i] - temp_R;
				OUT_I[i_lower] = X_I[i] - temp_I;
				OUT_R[i] = X_R[i] + temp_R;
				OUT_I[i] = X_I[i] + temp_I;
			}
			k+=step;
		}
		step=step/2;

}
//=======================END: FFT fft_stage_last =========================

