#include <stdio.h>
#include <stdlib.h>
#include "gen_laplace_mat.h"

int main(int argc, char *argv[]) {
	int n, i;
	int *diag, *l_upper, *l_lower, *u_upper, *u_lower;
	
	if (argc != 2) {
		printf("Need 1 arg\n");
		return 1;
	}
	
	n = atoi(argv[1]);
	
	gen_laplace_mat(n, &diag, &l_upper, &l_lower, &u_upper, &u_lower);
	
	printf("Diagonal:\n");
	for (i = 0; i < n * n; i++) {
		printf("%d ", diag[i]);
	}
	
	printf("\n\nU mat, lower band:\n");
	for (i = 0; i < n * n - 1; i++) {
		printf("%d ", u_lower[i]);
	}
	
	printf("\n\nL mat, upper band:\n");
	for (i = 0; i < n * n - 1; i++) {
		printf("%d ", l_upper[i]);
	}
	
	printf("\n\nU mat, upper band:\n");
	for (i = 0; i < n * n - n; i++) {
		printf("%d ", u_upper[i]);
	}
	
	printf("\n\nL mat, lower band:\n");
	for (i = 0; i < n * n - n; i++) {
		printf("%d ", l_lower[i]);
	}
	
	printf("\n");
	
	return 0;
}
