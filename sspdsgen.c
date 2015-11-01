#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>

void init_rand() {
  srand((unsigned)time(NULL));
}

double rand_within_range(double a, double b) {
  return a + (rand() / ( RAND_MAX / (b-a) ) ) ;
}

double ** random_sspdm(int n, double r, int *nz) {
  double **sspdm, tmp;
  int i, j;

  sspdm = (double **) malloc(n*sizeof(double*));
  assert(sspdm);
  for (i = 0; i < n; i++) { /* Creating a matrix of zeros. */
    sspdm[i] = (double *) calloc(n, sizeof(double));
    assert(sspdm[i]);
  }

  /* First we put 1 at each diagonal position. */
  for (i = 0; i < n; i++) {
    sspdm[i][i] = 1;
  }

  /* Now we put a random number from the uniform distribution on [-1, 1)
   * at each off-diagonal position (maintaining the symmetry A=A^T). 
   * Replacing each off-diagonal entry with a_ij if |a_ij| < r.
   */
  *nz = 0;
  for (i = 1; i < n; i++) {
    for (j = 0; j < i; j++) {
      tmp = rand_within_range(-1, 1);
      if (!(fabs(tmp) > r)) {
	(*nz)++;
	sspdm[i][j] = tmp;
	sspdm[j][i] = sspdm[i][j];
      }
    }
  }

  /* [DEBUG] Printing the matrix.
  printf("SSPDM:\n");
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      printf(" %f", sspdm[i][j]);
    }
    putchar('\n');
  }
  */

  return sspdm;
}

double * random_b(int n) {
  double *b;
  int i;

  b = (double *) calloc(n, sizeof(double));
  assert(b);
  for (i = 0; i < n; i++) {
    b[i] = rand_within_range(-10, 10);
  }

  /* [DEBUG] Printing vector b.
  puts("B:");
  for (i = 0; i < n; i++) {
  printf(" %f", i, b[i]);
  }
  putchar('\n');
  */
  
  return b;
}

int main(int argc, char **argv) {
  double **sspdm, *b, r;
  int i, j, n, nz;

  if (argc < 2) {
    printf("Few parameters! Please inform (at least) the dimension of the system.\n");
  } else if ((argc >= 2) && (argc <= 3)) {
    if (argc == 2) {
      n = atoi(argv[argc-1]);
      r = 0.01;
    } else {
      n = atoi(argv[argc-2]);
      r = atof(argv[argc-1]);
    }
    /* Generating the sparse system. */
    init_rand();  
    sspdm = random_sspdm(n, r, &nz);
    b = random_b(n);
    
    /* Printing to stdout the generated sparse system. */
    printf("%d\n", nz*2 + n);
    printf("%d\n", n);
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	printf("%d %d %f\n", i, j, sspdm[i][j]);
      }
    }
    for (i = 0; i < n; i++) {
      printf("%d %f\n", i, b[i]);
    }
    
    /* Releasing memory. */
    for(i = 0; i < n; i++) {
      free(sspdm[i]);
    }
    free(sspdm);
    free(b);
  } else {
    printf("Too many parameters!\n");
  }

  return 0;
}
