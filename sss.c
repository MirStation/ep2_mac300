#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

#include "sss.h"

void init_rand() {
  srand((unsigned)time(NULL));
}

double rand_within_range(double a, double b) {
  return a + (rand() / ( RAND_MAX / (b-a) ) ) ;
}

void random_spdm(int n, double r) {
  double **spdm, tmp;
  int i, j, k;

  spdm = (double **) malloc(n*sizeof(double*));
  assert(spdm);
  for (i = 0; i < n; i++) { /* Creating a matrix of zeros. */
    spdm[i] = (double *) calloc(n, sizeof(double));
    assert(spdm[i]);
  }

  /* First we put 1 at each diagonal position. */
  for (i = 0; i < n; i++) {
    spdm[i][i] = 1;
  }

  /* Now we put a random number from the uniform distribution on [-1, 1]
   * at each off-diagonal position (maintaining the symmetry A=A^T). 
   * Replacing each off-diagonal entry with a_ij if |a_ij| < r.
   */
  init_rand();
  k = 0;
  for (i = 1; i < n; i++) {
    for (j = 0; j < i; j++) {
      tmp = rand_within_range(-1, 1);
      if (!(fabs(tmp) > r)) {
	k++;
	printf("tmp: %f\n", tmp);
	spdm[i][j] = tmp;
	spdm[j][i] = spdm[i][j];
      }
    }
  }

  /* Printing the matrix. */
  printf("SPDM:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      printf(" %f", spdm[i][j]);
    }
    putchar('\n');
  }

  printf("nnz: %d\n", (k*2) + n);
  
  /* Releasing memory. */
  for(i = 0; i < n; i++) {
    free(spdm[i]);
  }
  free(spdm);
}

void sss_spm_v(SSS A, double *x, double **y) {
  int r, c, j;
  for (r = 0; r < A->n; r++) {
    (*y)[r] = A->dvalues[r] * x[r];
    for (j = A->rowptr[r]; j < A->rowptr[r + 1]; r++) {
      c = A->colind[j];
      (*y)[r] = (*y)[r] + A->values[j] * x[c];
      (*y)[c] = (*y)[c] + A->values[j] * x[r];
    }
  }
}
