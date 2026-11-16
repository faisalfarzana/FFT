
#include"fft_stages.h"
#include "math.h"
/*
 * Function: fft_stage_first_op calculates the first stage of FFT.
 * 
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal
	stage: stage number: //Possible values of STAGE=2,3,4,5,6,7,8,9

OUTPUT:
	OUT_R, OUT_I[]: Real and Imag parts of Complex signal
 */

//=======================BEGIN: FFT fft_stages =========================
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int stage, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

	//Write a code that computes any arbitary stages of the FFT
	DTYPE temp_R;		/*temporary storage complex variable*/
	DTYPE temp_I;		/*temporary storage complex variable*/


	int i,j,k;			/* loop indexes */
	int i_lower;		/* Index of lower point in butterfly */
	int step;

	int DFTpts;
	int numBF;			/*Butterfly Width*/

	int N2 = SIZE2;	/* N2=N>>1 */


	/*=======================BEGIN: FFT=========================*/
	// Do M stages of butterflies
	step=N2;
	DTYPE a, e, c, s;
	DFTpts = 1 << stage;		// DFT = 2^stage = points in sub DFT
	numBF = DFTpts/2; 			// Butterfly WIDTHS in sub-DFT
	k=0;

	e = -6.283185307178/DFTpts;

	a = 0.0;
		// Perform butterflies for j-th stage
	butterfly:for(j=0; j<numBF; j++)
		{

			c = cos(a);
			s = sin(a);
			a = a + e;

			// Compute butterflies that use same W**k
			DFTpts:for(i=j; i<SIZE; i += DFTpts)
			{

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
//=======================END: FFT fft_stages =========================

