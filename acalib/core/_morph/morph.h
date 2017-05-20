#ifndef MORPHOLOGY_HEADER
#define MORPHOLOGY_HEADER

void differenceImpl(double* cumulativeSum, double* difference, int n);
void segmentationImpl(double* diff, double* boxing, int n);
void erosionImpl(double* boxing, double* blocking, int n);

#endif
