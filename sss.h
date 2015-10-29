#ifndef SPARSE_SYMMETRIC_SKYLINE_H
#define SPARSE_SYMMETRIC_SKYLINE_H

struct SparceSymmetricSkyline {
  double *dvalues;
  int *rowptr;
  int *colind;
  double *values;
  int n;
  int nnz;
};
typedef struct SparceSymmetricSkyline * SSS;

void random_spdm(int n, double r);

void sss_spm_v(SSS A, double *x, double **y);

#endif
