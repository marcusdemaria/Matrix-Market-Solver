#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "functions.h"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Handle the inputs
    const char *filename = argv[1];
    CSRMatrix A;
    ReadMMtoCSR(filename, &A);

    // Initializing all the vector b (in Ax=b)
    double *b = (double *)malloc(A.num_rows * sizeof(double));
    
    // Set all elements of b to 1
    for (int i = 0; i < A.num_rows; ++i) {
        b[i] = 1.0;
    }

    // Timing the solver
    clock_t start_time = clock();

    // Call your solver function
    double *x = (double *)malloc(A.num_rows * sizeof(double)); //this will store the solution
    int max_iterations = 1000000; //pow(A.num_rows,4);  //adjust based on size of matrix
    double tolerance = 1e-14;  //adjust based on what convergence tolerance is

    // Call the solver
    solver(&A, b, x, max_iterations, tolerance);

    // Timing the solver
    clock_t end_time = clock();
    double cpu_time_taken = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;

    // Print out the matrix information
    printf("The matrix name: %s\n", filename);
    printf("The dimension of the matrix: %d by %d\n", A.num_rows, A.num_cols);
    printf("Number of non-zeros: %d\n", A.num_non_zeros);

    // Print out the CPU time taken to solve Ax = b
    printf("CPU time taken solve Ax=b: %f\n", cpu_time_taken);

    // Compute and print the residual norm
    double *residual = (double *)malloc(A.num_rows * sizeof(double));
    compute_residual(&A, b, x, residual);
    double residual_norm = compute_norm(residual, A.num_rows);

    printf("Residual Norm: %.16e\n", residual_norm);

    /*
    // Report the results in a tex file
    FILE *readme = fopen("README.tex", "w");
    if (readme != NULL) {
        fprintf(readme, "Results for %s:\n", filename);
        fprintf(readme, "The matrix dimensions: %d by %d\n", A.num_rows, A.num_cols);
        fprintf(readme, "Number of non-zeros: %d\n", A.num_non_zeros);
        fprintf(readme, "CPU time taken to solve Ax=b: %f\n", cpu_time_taken);
        fprintf(readme, "Residual Norm: %f\n", residual_norm);
        fclose(readme);
    } else {
        fprintf(stderr, "Failed to open README.tex for writing.\n");
    }
    */

    // Free allocated memory
    free(b);
    free(x);
    free(residual);

    free(A.col_ind);
    free(A.row_ptr);
    free(A.csr_data);

    return EXIT_SUCCESS;
}
