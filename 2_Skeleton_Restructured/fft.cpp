/*
This is traditional 2-radix DIT FFT algorithm implementation.
INPUT:
	In_R, In_I[]: Real and Imag parts of Complex signal

OUTPUT:
	Out_R, Out_I[]: Real and Imag parts of Complex signal
*/

#include <math.h>

#include "fft.h"
#include <string.h>

void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE]);
void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);
void fft_stage_first(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int STAGES, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);
void fft_stage_last(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]);


void write_output(DTYPE Stage10_R[SIZE], DTYPE Stage10_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {
    for (int i = 0; i < SIZE; i++) {
        #pragma HLS PIPELINE II=1
        OUT_R[i] = Stage10_R[i];
        OUT_I[i] = Stage10_I[i];
    }
}

void fft(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE])
{

    #pragma HLS INTERFACE m_axi depth=1024 port=X_R offset=slave bundle=data0
    #pragma HLS INTERFACE m_axi depth=1024 port=X_I offset=slave bundle=data1
    #pragma HLS INTERFACE m_axi depth=1024 port=OUT_R offset=slave bundle=data2
    #pragma HLS INTERFACE m_axi depth=1024 port=OUT_I offset=slave bundle=data3

    #pragma HLS INTERFACE s_axilite port=X_R bundle=control
    #pragma HLS INTERFACE s_axilite port=X_I bundle=control
    #pragma HLS INTERFACE s_axilite port=OUT_R bundle=control
    #pragma HLS INTERFACE s_axilite port=OUT_I bundle=control
    #pragma HLS INTERFACE s_axilite port=return bundle=control

	//task-level pipekining across all submoduls
    #pragma HLS dataflow

	DTYPE buff_I[SIZE];
	DTYPE buff_R[SIZE];


	memcpy(buff_I, X_I, SIZE * sizeof(DTYPE));
	memcpy(buff_R, X_R, SIZE * sizeof(DTYPE));


	//Call fft
	DTYPE Stage0_R[SIZE], Stage0_I[SIZE];
	DTYPE Stage1_R[SIZE], Stage1_I[SIZE];
	DTYPE Stage2_R[SIZE], Stage2_I[SIZE];
	DTYPE Stage3_R[SIZE], Stage3_I[SIZE];
	DTYPE Stage4_R[SIZE], Stage4_I[SIZE];
	DTYPE Stage5_R[SIZE], Stage5_I[SIZE];
	DTYPE Stage6_R[SIZE], Stage6_I[SIZE];
	DTYPE Stage7_R[SIZE], Stage7_I[SIZE];
	DTYPE Stage8_R[SIZE], Stage8_I[SIZE];
	DTYPE Stage9_R[SIZE], Stage9_I[SIZE];
    DTYPE Stage10_R[SIZE], Stage10_I[SIZE];

	bit_reverse(buff_R, buff_I, Stage0_R, Stage0_I);

	fft_stage_first(Stage0_R, Stage0_I, Stage1_R, Stage1_I);
	fft_stages(Stage1_R, Stage1_I, 2, Stage2_R, Stage2_I);
	fft_stages(Stage2_R, Stage2_I, 3, Stage3_R, Stage3_I);
	fft_stages(Stage3_R, Stage3_I, 4, Stage4_R, Stage4_I);
	fft_stages(Stage4_R, Stage4_I, 5, Stage5_R, Stage5_I);
	fft_stages(Stage5_R, Stage5_I, 6, Stage6_R, Stage6_I);
	fft_stages(Stage6_R, Stage6_I, 7, Stage7_R, Stage7_I);
	fft_stages(Stage7_R, Stage7_I, 8, Stage8_R, Stage8_I);
	fft_stages(Stage8_R, Stage8_I, 9, Stage9_R, Stage9_I);
	
	

	DTYPE result_R[SIZE];
	DTYPE result_I[SIZE];

	fft_stage_last(Stage9_R, Stage9_I, Stage10_R, Stage10_I);
	write_output(Stage10_R, Stage10_I, result_R, result_I);
	


	memcpy(OUT_R, result_R, SIZE * sizeof(DTYPE));
	memcpy(OUT_I, result_I, SIZE * sizeof(DTYPE));

}

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
	// //Insert your code here
	// unsigned int reversed;
	// 	DTYPE temp_R, temp_I;

	// for(int x = 0; x < SIZE; x++)
	// {
	// 	#pragma HLS PIPELINE II=2
    //     // inter false tells HLS not to stall the pipeline worrying about memory overlaps
    //     #pragma HLS dependence variable=X_R type=inter false
	// 	//HLS warning flag clear for swpaping indices 
	// 	#pragma HLS dependence variable=X_I type=inter false
		
	// 	reversed = reverse_bits(x);

	// 	if(x < reversed)
	// 	{
	// 		// swap ral
    //         temp_R = X_R[x];
    //         X_R[x] = X_R[reversed];
    //         X_R[reversed] = temp_R;

    //         //swap Imaginary
    //         temp_I = X_I[x];
    //         X_I[x] = X_I[reversed];
    //         X_I[reversed] = temp_I;
	// 	}
	// }
}

//4 args instead of two, to handle dual
void bit_reverse(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]){
    for(int x = 0; x < SIZE; x++)
    {
        #pragma HLS PIPELINE II=1
        int reversed = reverse_bits(x);
        
        //1 read, 1 write per cycle)
        OUT_R[reversed] = X_R[x];
        OUT_I[reversed] = X_I[x];
    }
}


/*=======================BEGIN: FFT=========================*/
//stage 1
void fft_stage_first(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

	//Insert your code here
	DTYPE temp_R;		/*temporary storage complex variable*/
	DTYPE temp_I;		/*temporary storage complex variable*/


	int i;
	int i_lower;		/* Index of lower point in butterfly */
	//int step;
	int DFTpts = 2;
    int numBF = 1;

    DFTpts_loop:for(i = 0; i < SIZE; i += DFTpts)
    {
        #pragma HLS PIPELINE II=1
        
        i_lower = i + numBF; // i_lower is always i + 1
        
        //twiddle factor is 1 for the first stage
        temp_R = X_R[i_lower]; 
        temp_I = X_I[i_lower];

        OUT_R[i_lower] = X_R[i] - temp_R;
        OUT_I[i_lower] = X_I[i] - temp_I;
        OUT_R[i] = X_R[i] + temp_R;
        OUT_I[i] = X_I[i] + temp_I;
    }
}

//stages
void fft_stages(DTYPE X_R[SIZE], DTYPE X_I[SIZE], int stage, DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {

    int DFTpts = 1 << stage;
    int numBF = DFTpts / 2;
    int step = SIZE >> stage; //map w imag
    
    // Flatten into a single loop to prevent pipeline flushing overhead.
    // SIZE / 2 is always exactly 512 for a 1024-point FFT.
    butterfly_loop:for(int iter = 0; iter < (SIZE / 2); iter++) {
        #pragma HLS PIPELINE II=1
        // Tell Vivado the calculated indices will never overlap to guarantee II=2
        #pragma HLS dependence variable=OUT_R type=inter false
        #pragma HLS dependence variable=OUT_I type=inter false
        
        // Because 'stage' is inlined, numBF is a constant power of 2.
        // Vivado will optimize this modulo and division into zero-cycle bit-shifts.
        int j = iter & (numBF - 1);
        int group = iter >> (stage-1);
        
        int i = (group << stage) + j;
        int i_lower = i + numBF;
        
        //precomputed twiddle factor from fft.h LUT instead of cos()/sin()
        int twiddle_idx = j * step;
        DTYPE c = W_real[twiddle_idx];
        DTYPE s = W_imag[twiddle_idx];
        
        DTYPE temp_R = X_R[i_lower]*c - X_I[i_lower]*s;
        DTYPE temp_I = X_I[i_lower]*c + X_R[i_lower]*s;
        
        OUT_R[i_lower] = X_R[i] - temp_R;
        OUT_I[i_lower] = X_I[i] - temp_I;
        OUT_R[i] = X_R[i] + temp_R;
        OUT_I[i] = X_I[i] + temp_I;
    }

}


//last stage
void fft_stage_last(DTYPE X_R[SIZE], DTYPE X_I[SIZE], DTYPE OUT_R[SIZE], DTYPE OUT_I[SIZE]) {
    DTYPE temp_R;		/*temporary storage complex variable*/
	DTYPE temp_I;		/*temporary storage complex variable*/
   
    //force amath to LUT fabric to save ~40 DSPs and fit under the 220 limit
    #pragma HLS bind_op variable=temp_R op=fmul impl=fabric
    #pragma HLS bind_op variable=temp_I op=fmul impl=fabric
    #pragma HLS bind_op variable=temp_R op=fadd impl=fabric
    #pragma HLS bind_op variable=temp_I op=fadd impl=fabric

//Write code that computes last stage of FFT
	//Insert your code here
	


	int i,j;			/* loop indexes */
    int i_lower;		/* Index of lower point in butterfly */

	//last stage (1024 points): numBF = 512
    int numBF = SIZE / 2; 

    butterfly:for(j = 0; j < numBF; j++)
    {
        #pragma HLS PIPELINE II=1
        
        //inner loop is gone. i is just equal to j in the last stage.
        i_lower = j + numBF; 
        
        temp_R = X_R[i_lower]*W_real[j] - X_I[i_lower]*W_imag[j];
        temp_I = X_I[i_lower]*W_real[j] + X_R[i_lower]*W_imag[j];

        OUT_R[i_lower] = X_R[j] - temp_R;
        OUT_I[i_lower] = X_I[j] - temp_I;
        OUT_R[j] = X_R[j] + temp_R;
        OUT_I[j] = X_I[j] + temp_I;
    }

}
/*=======================END: FFT=========================*/

