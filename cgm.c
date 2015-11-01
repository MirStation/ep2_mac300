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

void sss_spm_v(SSS A, double *x, double *y) {
  int r, c, j;
  for (r = 0; r < A->n; r++) {
    y[r] = A->dvalues[r] * x[r];
    for (j = A->rowptr[r]; j < A->rowptr[r + 1]; j++) {
      c = A->colind[j];
      y[r] = y[r] + A->values[j] * x[c];
      y[c] = y[c] + A->values[j] * x[r];
    }
  }
}

SSS new_sss(double **A, int n, int nz) {
  SSS new;
  int i, j, k, l, first_row_nz_value;

  new = (SSS) malloc(sizeof(struct SparceSymmetricSkyline));
  assert(new);
  new->dvalues = (double *) malloc(n * sizeof(double));
  assert(new->dvalues);
  new->rowptr = (int *) malloc((n + 1) * sizeof(int));
  assert(new->rowptr);
  for (i = 0; i < (n + 1); i++) {
    new->rowptr[i] = -1;
  }
  new->colind = (int *) malloc(((nz - n) / 2) * sizeof(int));
  assert(new->colind);
  new->values = (double *) malloc(((nz - n) / 2) * sizeof(double));
  assert(new->values);
  new->n = n;
  new->nz = nz;
  
  k = 0;
  for (i = 0; i < n; i++) {
    first_row_nz_value = 1;
    for (j = 0; j < i; j++) {
      if (A[i][j] != 0) {
	new->colind[k] = j;
	new->values[k] = A[i][j];
	if (first_row_nz_value == 1) {
	  first_row_nz_value = 0;
	  new->rowptr[i] = k;
	  l = i - 1;
	  while((new->rowptr[l] == -1) && (l >= 0)) {
	    new->rowptr[l--] = new->rowptr[i];
	  }
	}
	k++;
      }
    }
    new->dvalues[i] = A[i][j];
  }
  new->rowptr[i] = (nz - n) / 2;

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

void vvmul(int n, double *v1, double *v2, double *res) {
  int i;
  (*res) = 0;
  for (i = 0; i < n; i++) {
    (*res) += v1[i] * v2[i];
  }
}

void cvmul(int n, double c, double *v, double *res) {
  int i;
  for (i = 0; i < n; i++) {
    res[i] = c * v[i];
  }
}

void vvsum(int n, double *v1, double *v2, double *res) {
  int i;
  for (i = 0; i < n; i++) {
    res[i] = v1[i] + v2[i];
  }
}

void vvsub(int n, double *v1, double *v2, double *res) {
  int i;
  for (i = 0; i < n; i++) {
    res[i] = v1[i] - v2[i];
  }
}

void cg(SSS A, double *b, double e, double *x) {
  double *p, *Ap, alfa, beta, delta_0, delta, delta_old, aux, *aux_v;
  int i, j;
  
  Ap = (double *) calloc(A->n, sizeof(double));
  assert(Ap);

  p = (double *) malloc(A->n * sizeof(double));
  assert(p);
  for (i = 0; i < A->n; i++) {
    p[i] = b[i];
  }

  aux_v = (double *) calloc(A->n, sizeof(double));
  assert(aux_v);

  vvmul(A->n, b, b, &delta);
  delta_0 = delta;

  i = 0;
  while ((i < A->n) && (delta > (pow(e,2) * delta_0))) {
    sss_spm_v(A, p, Ap);
    /*puts("Ap:");
    for (j = 0; j < A->n; j++) {
      printf(" %f\n", Ap[j]);
      }
    */
    vvmul(A->n, p, Ap, &aux);
    /*printf("aux: %f\n",aux);*/
    alfa =  delta / aux;
    /*printf("alfa: %f\n",alfa);*/

    cvmul(A->n, alfa, p, aux_v);
    /*puts("x-aux_v:");
    for (j = 0; j < A->n; j++) {
      printf(" %f\n", aux_v[j]);
      }*/
    vvsum(A->n, x, aux_v, x);
    /*puts("x:");
    for (j = 0; j < A->n; j++) {
      printf(" %f\n", x[j]);
      }*/
    
    cvmul(A->n, alfa, Ap, aux_v);
    /*puts("b-aux_v:");
    for (j = 0; j < A->n; j++) {
      printf(" %f\n", aux_v[j]);
      }*/
    vvsub(A->n, b, aux_v, b);
    /*puts("b:");
    for (j = 0; j < A->n; j++) {
      printf(" %f\n", b[j]);
      }*/

    delta_old = delta;
    vvmul(A->n, b, b, &delta);
    /*printf("delta: %f\n",delta);
      printf("delta_old: %f\n",delta_old);*/
    beta = delta / delta_old;
    /*printf("beta: %f\n",beta);*/
    
    cvmul(A->n, beta, p, aux_v);
    /*puts("p-aux_v:");
    for (j = 0; j < A->n; j++) {
      printf(" %f\n", aux_v[j]);
      }*/
    vvsum(A->n, b, aux_v, p);
    /*puts("p:");
    for (j = 0; j < A->n; j++) {
      printf(" %f\n", p[j]);
      }*/
    i++;
  }

  /*
  if (i < A->n) {
    printf("x^T = [");
    for (i = 0; i < n; i++) {
      printf(" %f", x[i]);
    }
    printf("]\n");
  } else {
    printf("Conjugate gradient method failed!");
  }
  */
  
  free(Ap);
  free(p);
  free(aux_v);
}

int main(int argc, char **argv) {
  SSS sss;
  double **A, *b, *x, v;
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

  x = (double *) calloc(n, sizeof(double));
  assert(x);
  
  cg(sss, b, 0.0001, x);

  puts("x^T:");
  for (i = 0; i < n; i++) {
    printf("%f\n", x[i]);
  }

  delete_sss(sss);
  for (i = 0; i < n; i++) {
    free(A[i]);
  }
  free(A);
  free(b);
    
  return 0;
}

