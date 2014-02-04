#include <stdlib.h>
#include <omp.h>
#include "gen_laplace_mat.h"

void stencil(int, int, int *, int *);

// Generates a 2D laplacian matrix that is
// n^2 x n^2
// O(N) algorithm
void gen_laplace_mat(int n, int **diag, int **l_upper_band, int **l_lower_band, int **u_upper_band, int **u_lower_band) {
	int diag_len = n * n;
	int mid_band_len = diag_len - 1;
	int outer_band_len = diag_len - n;
	int i, low_fill, up_fill;
	
	*diag = malloc(diag_len * sizeof(int));
	*l_upper_band = malloc(mid_band_len * sizeof(int));
	*u_lower_band = malloc(mid_band_len * sizeof(int));
	*l_lower_band = malloc(outer_band_len * sizeof(int));
	*u_upper_band = malloc(outer_band_len * sizeof(int));
	
#pragma omp parallel for private(low_fill, up_fill)
	for (i = 0; i < diag_len; i++) {
		(*diag)[i] = -4;
		
		stencil(n, i + 1, &low_fill, &up_fill);
		if (i > 0)
			(*l_upper_band)[i - 1] = low_fill;
		if (i < mid_band_len)
			(*u_lower_band)[i] = up_fill;
		
		if (i < outer_band_len) {
			(*l_lower_band)[i] = 1;
			(*u_upper_band)[i] = 1;
		}
	}
	
	return;
}

void stencil(int n, int center_pos, int *lower_fill, int *upper_fill) {
	*lower_fill = center_pos % n == 1 ? 0 : 1;
	*upper_fill = center_pos % n == 0 ? 0 : 1;
		
	return;
}
