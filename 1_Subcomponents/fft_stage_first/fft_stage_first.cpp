#include"fft_stage_first.h"
/*
 * Function: fft_stage_first_op calculates the first stage of the FFT.
 * 
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal

OUTPUT:
	OUT_R, OUT_I[]: Real and Imag parts of Complex signal
 */

//=======================BEGIN: FFT fft_stage_first_op =========================
void fft_stage_first(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

	//Insert your code here
	DTYPE temp_R;		/*temporary storage complex variable*/
	DTYPE temp_I;		/*temporary storage complex variable*/


	int i,j,k;			/* loop indexes */
	int i_lower;		/* Index of lower point in butterfly */
	int step;

	int stage = 1;
	int DFTpts;
	int numBF;			/*Butterfly Width*/

	int N2 = SIZE2;	/* N2=N>>1 */
	step=N2;
	DTYPE e;

		DFTpts = 1 << stage;		// DFT = 2^stage = points in sub DFT
		numBF = DFTpts/2; 			// Butterfly WIDTHS in sub-DFT
		k=0;

		e = -6.283185307178/DFTpts;
		// Perform butterflies for j-th stage
		butterfly:for(j=0; j<numBF; j++)
		{

			// Compute butterflies that use same W**k
			DFTpts:for(i=j; i<SIZE; i += DFTpts)
			{

				i_lower = i + numBF;			//index of lower point in butterfly
				temp_R = X_R[i_lower]*1;
				temp_I = X_I[i_lower]*1;

				OUT_R[i_lower] = X_R[i] - temp_R;
				OUT_I[i_lower] = X_I[i] - temp_I;
				OUT_R[i] = X_R[i] + temp_R;
				OUT_I[i] = X_I[i] + temp_I;
			}
			k+=step;
		}
		step=step/2;
	}

//=======================END: FFT fft_stage_first_op =========================

