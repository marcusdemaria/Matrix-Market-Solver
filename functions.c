#include "functions.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h> // Include for DBL_EPSILON
#include <string.h> // Include for symmetry

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix)
{
    // Open the file
    FILE *file = fopen(filename, "r"); //in read mode
    char line[1000]; //saving bits/memory for the file

    // Check for file opening error
    if (file == NULL)
    {
        printf("Error: Unable to open file %s\n", filename);
        return 1;
    }

    // Check if the matrix is symmetric
    int symmetric;
    fgets(line, sizeof(line), file);
    symmetric = (strstr(line, "symmetric") != NULL) ? 1 : 0;
    //if string "symmetric" found in first line, returns pointer to first occurance
    //ternary operate evaluates to true if symmetric is found (aka is not NULL) or else false

    // Skip comments and blank lines
    do
    {
        fgets(line, sizeof(line), file); //reads line from file and stores in array line, fgets reading up to size of line
    } while (line[0] == '%' || line[0] == '\n'); //checks if first character of each line is a % or newline character to skip

    // Scan the first line to get values needed
    sscanf(line, "%d %d %d", &(matrix->num_rows), &(matrix->num_cols), &(matrix->num_non_zeros));

    // Adjust non-zero count for symmetric matrix
    if (symmetric)
    {
        matrix->num_non_zeros += (matrix->num_non_zeros - matrix->num_cols);
    }

    // Dynamically allocate memory for matrix values
    int *temp_row = (int *)malloc((matrix->num_non_zeros) * sizeof(int));
    matrix->col_ind = (int *)malloc(matrix->num_non_zeros * sizeof(int));
    matrix->csr_data = (double *)malloc(matrix->num_non_zeros * sizeof(double));

    // Check for allocation failure
    if (temp_row == NULL || matrix->col_ind == NULL || matrix->csr_data == NULL)
    {
        printf("Error: Memory allocation failed.\n");
        return 0;
    }

    int row, col;
    double val;

    // Assign values from the file in CSR format
    for (int i = 0; i < matrix->num_non_zeros; ++i) //assigns values in file to correct array
    {
        fscanf(file, "%d %d %lf\n", &row, &col, &val);

        matrix->csr_data[i] = val;
        matrix->col_ind[i] = col - 1;

        temp_row[i] = row - 1;

        // Adjust for symmetric matrix
        if (symmetric && (matrix->col_ind[i] != temp_row[i]))
        {
            temp_row[i + 1] = matrix->col_ind[i];
            matrix->col_ind[i + 1] = temp_row[i];
            matrix->csr_data[i + 1] = matrix->csr_data[i];
            i++;
        }
    }
    fclose(file);

    int i, j, temp, min_index;
    double temp_csr;

    // Sorting algorithm
    // sorts based on rows from smallest to greatest
    for (i = 0; i < (matrix->num_non_zeros - 1); i++)
    {
        min_index = i; //assume current first value is smallest
        for (j = i + 1; j < (matrix->num_non_zeros); j++) //searching through each value of array to sort array
        {
            if (temp_row[j] < temp_row[min_index]) 
            {
                min_index = j; //if current row value less than assigned one, make it the new assigned one
            }
        }

        temp = temp_row[min_index]; //temporarily stores value in variable temp
        temp_row[min_index] = temp_row[i]; //replaces min index value with value at index i
        temp_row[i] = temp; //puts temp value into front of temp row, so i in temp row is the smallest row val

        temp = matrix->col_ind[min_index];
        matrix->col_ind[min_index] = matrix->col_ind[i];
        matrix->col_ind[i] = temp;

        temp_csr = matrix->csr_data[min_index];
        matrix->csr_data[min_index] = matrix->csr_data[i];
        matrix->csr_data[i] = temp_csr;
    }

    int ind = 1;
    int val_count = 1;

    // Convert rows array to the proper array
    for (int i = 1; i < matrix->num_non_zeros; i++)
    {
        if (temp_row[i] != temp_row[i - 1]) //checks if the temp row row values are different
        {
            temp_row[ind] = i; //if different, a new row is starting so it updates temp row to store current index
            ind++; 
            val_count++;
            //increment to move to next position in temp row
        }
    }
    val_count++; //is updated every time a new row starts, but needs to update one more time to match index format

    temp_row[ind] = matrix->num_non_zeros; //sets to num non zeros to end the array (last index is total)

    // Reallocate the rowptr array
    matrix->row_ptr = (int *)malloc((val_count) * sizeof(int)); //allocate to starting indicies of each row in sorted matrix
    for (int i = 0; i < val_count; i++)
    {
        matrix->row_ptr[i] = temp_row[i]; //copy values from temp row to the matrix row ptr
    }

    matrix->num_non_zeros -= (matrix->num_non_zeros / 2 - matrix->num_cols / 2); //adjust num non zeros back to how many in file reading
    free(temp_row); //free temp row memory

    return 1;
}

void spmv_csr(const CSRMatrix *A, const double *x, double *y) {
    // Implement sparse matrix-vector multiplication for CSR format

    // Iterate over each row of the matrix
    for (int i = 0; i < A->num_rows; i++){
        y[i]=0.0; // initalizes each element of result vector y to zero, that stores result of function
    }
    for (int i = 0; i < A->num_rows; i++){ //iterate over each row
        int start_row = A->row_ptr[i]; //start and end indicies
        int end_row = A->row_ptr[i+1];
        for (int j = start_row; j < end_row; j++){  //iterate from start to end
            y[i] += A->csr_data[j] *x[A->col_ind[j]]; //product of matrix entries and element of input vector x into each y value
        }
    }
}

void solver(const CSRMatrix *A, const double *b, double *x, int max_iterations, double tolerance) {
    int n = A->num_rows;

    //solution vector with zeros
    for (int i = 0; i < n; i++){
        x[i] = 1;
    }

    //Temp vector to store updated solution
    double *tempv = (double *)malloc(n * sizeof(double));
    //start iterations
    int iterate = 0;
    double max_change = 0;
    do
    {
        max_change = 0;

        //solve for updated solution
        for (int rowIndex = 0; rowIndex < n; rowIndex++){ //keeps track of row currently being solved
            double updatedValue = 0;
            double diagonalValue; //temporary values for calculations
            int rowStart = A->row_ptr[rowIndex];
            int rowEnd = A->row_ptr[rowIndex + 1]; //start and end of the current row

            for (int entryIndex = rowStart; entryIndex < rowEnd; entryIndex++){ //iterate through row
                if (A->col_ind[entryIndex] == rowIndex){
                    diagonalValue = A->csr_data[entryIndex]; //if current value has same row and col, store in diagonal value
                }
                else{
                    updatedValue += A->csr_data[entryIndex] * x[A->col_ind[entryIndex]]; //else produce product of matrix value and value in soln vector
                }
            }

            // Compute the new value for the current element
            double newValue = (b[rowIndex] - updatedValue) / diagonalValue; //difference between b and sum of other matrix elements multiplied by soln vec values

            // Update the solution vector and track the maximum change
            max_change = fmax(max_change, fabs(newValue - x[rowIndex])); //check for change in current element is largest change seen so far

            // Update the solution vector immediately in the same iteration
            x[rowIndex] = newValue; //update solution vector with newly calculate value from Jacobi formula to use in next iteration
        }

        // Print the max change after each iteration
        /*
        printf("Change after iteration %d: %e\n", iterate + 1, max_change);
        */

        /*
        //updated solution copied into old one
        for (int i = 0; i < n; i++){
            x[i] = tempv[i];
        }
        */

        iterate++;
    } while (max_change > tolerance && iterate < max_iterations);
    
    // Print what x vals are
    /*
    for (int i = 0; i < n; i++){
            printf("%f \n", x[i]);
    }
    free(tempv);
    */
}

void compute_residual(const CSRMatrix *A, const double *b, const double *x, double *residual) {
    // Compute the residual: residual = Ax - b
    spmv_csr(A, x, residual); //matrix multiplication called to get A*x
    for (int i = 0; i < A->num_rows; ++i) {
        residual[i] -= b[i]; //for every val, residual is A*x - b
    }
}

double compute_norm(const double *vector, int size) {
    double norm = 0.0; //initialize norm

    for (int i = 0; i < size; ++i) {
        norm += vector[i] * vector[i]; //for size of matrix, norm is vec * vec at each index
    }

    return sqrt(norm);
}
