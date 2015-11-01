#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

struct SparceSymmetricSkyline {
  double *dvalues;
  int *rowptr;
  int *colind;
  double *values;
  int n;
  int nz;
};
typedef struct SparceSymmetricSkyline * SSS;

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

SSS new_sss(double **A, int n, int nz) {
  SSS new;
  int i, j, k, l, first_row_nz_value;

  new = (SSS) malloc(sizeof(struct SparceSymmetricSkyline));
  assert(new);
  new->dvalues = (double *) malloc(n*sizeof(double));
  assert(new->dvalues);
  new->rowptr = (int *) malloc((n+1)*sizeof(int));
  assert(new->rowptr);
  new->colind = (int *) malloc(((nz - n) / 2)*sizeof(int));
  assert(new->colind);
  new->values = (double *) malloc(((nz - n) / 2)*sizeof(double));
  assert(new->values);
  new->n = n;
  new->nz = nz;

  k = 0;
  l = 0;
  for (i = 0; i < n; i++) {
    first_row_nz_value = 1;
    for (j = 0; j < i; j++) {
      if (A[i][j] != 0) {
	if (first_row_nz_value == 1) {
	  first_row_nz_value = 0;
	  new->rowptr[k++] = j;
	}
	new->colind[l] = j;
	new->values[l] = A[i][j];
	l++;
      }
    }
    new->dvalues[i] = A[i][j];
    if (first_row_nz_value == 1) {
      new->rowptr[k++] = j;
    }
  }
  new->rowptr[k] = (nz - n) / 2;

  /* [DEBUG] Printing SSS. 
  puts("dvalues:");
  for (i = 0; i < n; i++) {
    printf(" %f", new->dvalues[i]);
  }
  putchar('\n');
  puts("rowptr:");
  for (i = 0; i < (n + 1); i++) {
    printf(" %d", new->rowptr[i]);
  }
  putchar('\n');
  puts("colind:");
  for (i = 0; i < ((nz - n) / 2); i++) {
    printf(" %d", new->colind[i]);
  }
  putchar('\n');
  puts("values:");
  for (i = 0; i < ((nz - n) / 2); i++) {
    printf(" %f", new->values[i]);
  }
  putchar('\n');
  printf("n:\n %d\n", new->n);
  printf("nz:\n %d\n", new->nz);
  */
  
  return new;
}

void delete_sss(SSS sm) {
  free(sm->dvalues);
  free(sm->rowptr);
  free(sm->colind);
  free(sm->values);
  free(sm);
}

/* Do not use it.
void cg(int imax, int n, double **A, double *b, double e){
  double *r, *x, *p, af, bt, *s0, s1;
  int i;
  x = (double *) calloc(n, sizeof(double));
  assert(x);
  r = (double *) malloc(n * sizeof(double));
  assert(r);
  p = (double *) malloc(n * sizeof(double));
  assert(p);
  for (i = 1; i < n; i++) {
    p[i] = b[i];
  }
  for (i = 1; i < imax; i++) {
    s0 = Mv(n, A, p); 
    s1 = vvmul(n, b, b);
    af = s1 / vvmul(n, p, s0);
    x = vvsum(n, x, vemul(n, p, af));
    r = vvsub(n, b, vemul(n, s0, af));
    free(s0);
    if (sqrt(vvmul(n, r, r)) < e) break;
    bt = vvmul(n, r, r) / vvmul(n, b, b);
    p = vvsum(n, r, vemul(n, p, bt));
    b = r;
  }
  if (i < imax) {
    printf("x^T = [");
    for (i = 0; i < n; i++) {
      printf(" %f", x[i]);
    }
    printf("]\n");
  } else {
    printf("Conjugate gradient method failed!");
  }
  free(r);
  free(p);
}
*/

int main(int argc, char **argv) {
  SSS sss;
  double **A, *b, v;
  int i, j, k, nz, n;

  scanf("%d", &nz);
  scanf("%d", &n);

  A = (double **) malloc(n * sizeof(double*));
  assert(A);
  for (i = 0; i < n; i++) {
    A[i] = (double *) malloc(n * sizeof(double));
    assert(A[i]);
  }
  for (i = 0; i < (n*n); i++) {
    scanf("%d %d %lf", &j, &k, &v);
    A[j][k] = v;
  }
  
  b = (double *) malloc(n * sizeof(double));
  assert(b);
  for (i = 0; i < n; i++) {
    scanf("%d %lf", &j, &v);
    b[j] = v;
  }

  sss = new_sss(A, n, nz);

  delete_sss(sss);
  for (i = 0; i < n; i++) {
    free(A[i]);
  }
  free(A);
  free(b);
    
  return 0;
}

