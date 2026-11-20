#include "fft.h"
#include <stdio.h>

int closestSymbol(DTYPE R, DTYPE I)
{
	if (R >= 0 && I >= 0) { return 0; }
	else if (R < 0 && I >= 0) { return 1; }
	else if (R < 0 && I < 0) { return 2; }
	else if (R >= 0 && I < 0) { return 3; }
}

void qpsk_decode(DTYPE R[SIZE], DTYPE I[SIZE], int D[SIZE]) {

	//Write your code here
	
	for (int i = 0; i < SIZE; i++)
	{
		D[i] = closestSymbol(R[i], I[i]);
	}

}

