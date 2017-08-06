#include <string.h>
#include "morph.h"

void differenceImpl(double* cumulativeSum, double* difference, int n)
{
    difference[0] = cumulativeSum[0];
    int i;
    for(i = 1; i < n; i++)
    {
        difference[i] = cumulativeSum[i] - difference[i-1];
    }
}

void segmentationImpl(double* diff, double* boxing, int n)
{
    int i;
	for(i = 0; i < n; i++)
    {
        boxing[i] = 0;
    }
    for(i = 1; i < n-1; i++)
    {
        boxing[i] = 1;
        int A = diff[i] < diff[i-1];
        int B = diff[i] < diff[i+1];
        int C = diff[i] > diff[i-1];
        int D = diff[i] > diff[i+1];
        if((A && B) || (C && D))
        {
            boxing[i] = 0;
        }
    }
}

void erosionImpl(double* boxing, double* blocking, int n)
{
    int i;
    for(i = 0; i < n; i++)
    {
        blocking[i] = 0;
    }
    for(i = 1; i < n-1; i++)
    {
        blocking[i] = boxing[i];
        int A = boxing[i-1] == 0;
        int B = boxing[i] == 1;
        int C = boxing[i+1] == 0;
        if(A && B && C)
        {
            blocking[i] = 0;
        }
    }
    memcpy(boxing, blocking, n*sizeof(double));
    for(i = 1; i < n-1; i++)
    {
        int A = blocking[i-1] == 0;
        int B = blocking[i] == 1;
        int C = blocking[i+1] == 0;
        if(A && B)
        {
            boxing[i-1] = 1;
        }
        if(B && C)
        {
            boxing[i+1] = 1;
        }
    }
}
