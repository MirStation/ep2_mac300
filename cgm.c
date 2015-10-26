#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

double * vvsum(int n, double *v0, double *v1) {
  int i;
  for (i = 0; i < n; i++) {
    v0[i] += v1[i]; 
  }
  return v0;
}

double * vvsub(int n, double *v0, double *v1) {
  int i;
  for (i = 0; i < n; i++) {
    v0[i] -= v1[i]; 
  }
  return v0;
}

double * vemul(int n, double *v, double e) {
  int i;
  for (i = 0; i < n; i++) {
    v[i] *= e; 
  }
  return v;
}

double vvmul(int n, double *v0, double *v1) {
  int i;
  double r = 0;
  for (i = 0; i < n; i++) {
    r += v0[i] * v1[i]; 
  }
  return r;
}

double * Mv(int n, double **M, double *v) {
  int i;
  double *vr = (double *) calloc(n, sizeof(double));
  assert(vr);
  for (i = 0; i < n; i++) {
    vr[i] = vvmul(n, M[i], v); 
  }
  return vr;
}

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


int main(int argc, char **argv) {
  printf("Conjugate Gradient Method!\n");
  return 0;
}
