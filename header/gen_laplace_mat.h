#ifndef _GEN_LAPLACE_MAT_H
#define _GEN_LAPLACE_MAT_H

// Generates a 2D laplacian matrix that is
// n^2 x n^2
// O(N) algorithm
void gen_laplace_mat(int n, int **diag, int **l_upper_band, int **l_lower_band, int **u_upper_band, int **u_lower_band);

#endif
