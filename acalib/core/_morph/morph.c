#include "morph.h"

void differenceImpl(double* cumulativeSum, double* difference, int n)
{
    difference[0] = cumulativeSum[0];
    for(int i = 1; i < n; i++)
    {
        difference[i] = cumulativeSum[i] - difference[i-1];
    }
}

void segmentationImpl(double* diff, double* boxing, int n)
{
    for(int i = 1; i < n-1; i++)
    {
        if((diff[i] < diff[i-1] && diff[i] < diff[i+1]) || (diff[i] > diff[i-1] && diff[i] > diff[i+1]))
        {
            boxing[i] = 0;
        }
    }
}

void erosionImpl(double* boxing, double* blocking, int n)
{
    for(int i = 0; i < n-1; i++)
    {
        blocking[i] = boxing[i];
        if(boxing[i-1] == 0 && boxing[i] == 1 && boxing[i+1] == 0)
        {
            blocking[i] = 0;
        }
    }
    for(int i = 0; i < n; i++)
    {
        boxing[i] = blocking[i];
    }
    for(int i = 1; i < n-1; i++)
    {
        if(blocking[i-1] == 0 && blocking[i] == 1)
        {
            boxing[i-1] = 1;
        }
        if(blocking[i] == 1 && blocking[i+1] == 0)
        {
            boxing[i+1] = 1;
        }
    }
}
