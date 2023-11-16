#include "sparse_matrix_formats.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <ostream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring> // for std::memcpy

#include <iostream>
#include <fstream> // Include this header for file stream operations
#include <string>
#include <cstdlib>
#include <cmath> // Include this header for mathematical operations
#include <string> // Include this header for mathematical operations
#include <stdbool.h>


extern "C" {
void wrapper_get_upper_csr_bin_c(
    const char* filename, 
    double** a, 
    int** ia, 
    int** p_ja, 
    int* nnz,
    int* dim, 
    bool zero_based,
    double** b_rhs, double** sol_x, int* num_rhs);
}


extern "C" {
void get_upper_csr_bin_c_args(int argc, char **argv, double** a, int** ia, int** p_ja, int* nnz, int* dim, bool get_b, bool get_x, double** b_rhs, double** sol_x, int* num_rhs);
}

extern "C" {
void get_upper_csr_bin_c_filename(const char *filename, double** a, int** ia, int** p_ja, int* nnz, int* dim, bool get_b, bool get_x, double** b_rhs, double** sol_x, int* num_rhs);
}


void printMatrixInformation_sympack(const double* values, 
                            int nnz,
                            int dim, 
                            const int* sympack_csc_i, 
                            const int* sympack_csc_p_j, 
                            const int* sympack_x_order) {
    printf("Matrix dimension: %d\n", dim);
    printf("Number of non-zero elements: %d\n\n", nnz);

    printf("sympack_x (non-zero values):\n");
    for (int i = 0; i < nnz; ++i) {
        printf("%f ", values[i]);
    }
    printf("\n\n");

    printf("CSC index array (sympack_csc_i):\n");
    for (int i = 0; i < nnz; ++i) {
        printf("%d ", sympack_csc_i[i]);
    }
    printf("\n\n");

    printf("CSC column pointer array (sympack_csc_p_j):\n");
    for (int i = 0; i < dim + 1; ++i) {
        printf("%d ", sympack_csc_p_j[i]);
    }
    printf("\n\n");

    //printf("Sorted order array (sympack_x_order):\n");
    //for (int i = 0; i < nnz; ++i) {
    //    printf("%d ", sympack_x_order[i]);
    //}
    //printf("\n\n");
}

void printMatrixInformation_pexsi(const double* values, 
                            int nnz,
                            int dim, 
                            const int* pexsi_csc_i, 
                            const int* pexsi_csc_p_j, 
                            const int* pexsi_x_order) {

    int full_nnz = 2*nnz - dim;
    printf("Matrix dimension: %d\n", dim);
    printf("Number of non-zero elements: %d\n\n", full_nnz);

    printf("Array a (non-zero values):\n");
    for (int i = 0; i < full_nnz; ++i) {
        printf("%f ", values[i]);
    }
    printf("\n\n");

    printf("CSC index array (pexsi_csc_i):\n");
    for (int i = 0; i < full_nnz; ++i) {
        printf("%d ", pexsi_csc_i[i]);
    }
    printf("\n\n");

    printf("CSC column pointer array (pexsi_csc_p_j):\n");
    for (int i = 0; i < dim + 1; ++i) {
        printf("%d ", pexsi_csc_p_j[i]);
    }
    printf("\n\n");

    //printf("Sorted order array (pexsi_x_order):\n");
    //for (int i = 0; i < full_nnz; ++i) {
    //    printf("%d ", pexsi_x_order[i]);
    //}
    //printf("\n");
}


void update_x_sympack(const double* a, int nnz, int dim, const int* sympack_x_order, double* sympack_x){
    for (int i = 0; i < nnz; ++i) {
        sympack_x[i] = a[sympack_x_order[i]];
    }
}

void update_x_pexsi(const double* a, int nnz, int dim, const int* pexsi_x_order, double* pexsi_x){
    int full_nnz = 2*nnz - dim;
    for (int i = 0; i < full_nnz; ++i) {
        pexsi_x[i] = a[pexsi_x_order[i]];
    }
}

void print_x_sympack(const double* sympack_x, int nnz, int dim) {
    printf("\nSympack x values:\n");
    for (int i = 0; i < nnz; ++i) {
        printf("%f ", sympack_x[i]);  // sympack_x[i] is already a double, no dereferencing needed
    }
    printf("\n");
}

void print_x_pexsi(const double* pexsi_x, int nnz, int dim) {
    int full_nnz = 2 * nnz - dim;
    printf("\nPexsi x values:\n");
    for (int i = 0; i < full_nnz; ++i) {
        printf("%f ", pexsi_x[i]);  // pexsi_x[i] is already a double, no dereferencing needed
    }
    printf("\n\n");
}


extern "C" {
    int inla_format_interface_sympack(double** a, int** ia, int** p_ja, int* nnz, int* dim,
                            int** sympack_csc_i, int** sympack_csc_p_j, int** sympack_x_order,
                            double** sympack_x, int* num_calls_sympack);
}

extern "C" {
    int inla_format_interface_pexsi(double** a, int** ia, int** p_ja, int* nnz, int* dim,
                            int** pexsi_csc_i, int** pexsi_csc_p_j, int** pexsi_x_order,
                            double** pexsi_x, int* num_calls_pexsi);
}

extern "C" {
    int inla_sympack_pexsi_format_interface(
        double* a, int* ia, int* p_ja, int nnz, int dim,
        int* sympack_csc_i, int* sympack_csc_p_j, int* sympack_x_order, double* sympack_x, 
        int* pexsi_csc_i, int* pexsi_csc_p_j, int* pexsi_x_order, double* pexsi_x, 
        int* num_calls_sympack, int* num_calls_pexsi, bool sympack_interface, bool pexsi_interface
    ) {

        
        if(sympack_interface) inla_format_interface_sympack(&a, &ia, &p_ja, &nnz, &dim, 
                            &sympack_csc_i, &sympack_csc_p_j, &sympack_x_order, &sympack_x,
                            num_calls_sympack);

        if(pexsi_interface) inla_format_interface_pexsi(&a, &ia, &p_ja, &nnz, &dim, 
                            &pexsi_csc_i, &pexsi_csc_p_j, &pexsi_x_order, &pexsi_x,
                            num_calls_pexsi);

        return 1; 
    }
}


int original_main(int argc, char **argv) {

    //----------------------------------------------------
    //               phase 1: input
    //----------------------------------------------------

    double* a = NULL;
    int* ia = NULL;
    int* p_ja = NULL;
    int nnz;
    int dim;
    int num_rhs;
    double *rhs_b, *sol_x;

    get_upper_csr_bin_c_args(argc, argv, &a, &ia, &p_ja, &nnz, &dim,false, false, &rhs_b, &sol_x, &num_rhs);

    //----------------------------------------------------
    //               phase 2: output
    //----------------------------------------------------

    int* sympack_csc_i = NULL;
    int* sympack_csc_p_j = NULL;
    int* sympack_x_order = NULL;
    int num_calls_sympack = 0;
    double* sympack_x = (double*)malloc(nnz * sizeof(double));


    int* pexsi_csc_i = NULL;
    int* pexsi_csc_p_j = NULL;
    int* pexsi_x_order = NULL;
    int num_calls_pexsi = 0;
    double* pexsi_x = (double*)malloc((2*nnz-dim) * sizeof(double));
    

    //first call
    inla_format_interface_sympack(&a, &ia, &p_ja, &nnz, &dim, 
                          &sympack_csc_i, &sympack_csc_p_j, &sympack_x_order, &sympack_x,
                          &num_calls_sympack);

    //first call
    inla_format_interface_pexsi(&a, &ia, &p_ja, &nnz, &dim, 
                        &pexsi_csc_i, &pexsi_csc_p_j, &pexsi_x_order, &pexsi_x,
                        &num_calls_pexsi);


    //second call
    inla_format_interface_sympack(&a, &ia, &p_ja, &nnz, &dim, 
                          &sympack_csc_i, &sympack_csc_p_j, &sympack_x_order, &sympack_x,
                          &num_calls_sympack);
    
    //second call
    inla_format_interface_pexsi(&a, &ia, &p_ja, &nnz, &dim, 
                        &pexsi_csc_i, &pexsi_csc_p_j, &pexsi_x_order, &pexsi_x,
                        &num_calls_pexsi);

    printMatrixInformation_sympack(sympack_x, nnz, dim, sympack_csc_i, sympack_csc_p_j, sympack_x_order);
    printMatrixInformation_pexsi(pexsi_x, nnz, dim, pexsi_csc_i, pexsi_csc_p_j, pexsi_x_order);


    // Free allocated memory before exiting
    free(pexsi_x);
    free(sympack_x);
    free(sympack_csc_i);
    free(sympack_csc_p_j);
    free(sympack_x_order);
    free(pexsi_csc_i);
    free(pexsi_csc_p_j);
    free(pexsi_x_order);
    free(a);
    free(ia);
    free(p_ja);

    return 0; // Return 0 to indicate success
}


int inla_format_interface_sympack(double** a, int** ia, int** p_ja, int* nnz, int* dim,
                          int** sympack_csc_i, int** sympack_csc_p_j, int** sympack_x_order,
                          double** sympack_x, int* num_calls_sympack) {

    //matrix.export_rb_format(std::string(argv[1]));

    if ((*num_calls_sympack) == 0) {
        SymmetricSparseMatrixFormats matrix(*a, *ia, *p_ja, *nnz, *dim);
        //matrix.print_upper_csr();

        matrix.get_lower_csc_from_upper_csr(true);
        matrix.exportToGlobalLowerCSCPointers(*sympack_csc_i, *sympack_csc_p_j, *sympack_x_order);
        (*num_calls_sympack)++;

    }else{

        update_x_sympack(*a, *nnz, *dim, *sympack_x_order, *sympack_x);
        //print_x_sympack(*sympack_x, *nnz, *dim); // Dereference sympack_x here
    }

    return 0; // Return 0 to indicate success

}

int inla_format_interface_pexsi(double** a, int** ia, int** p_ja, int* nnz, int* dim,
                          int** pexsi_csc_i, int** pexsi_csc_p_j, int** pexsi_x_order,
                          double** pexsi_x, int* num_calls_pexsi) {

    if ((*num_calls_pexsi) == 0) {
        SymmetricSparseMatrixFormats matrix(*a, *ia, *p_ja, *nnz, *dim);
        //matrix.print_upper_csr();

        matrix.get_full_csc_from_upper_csr(false);
        matrix.exportToGlobalFullCSCPointers(*pexsi_csc_i, *pexsi_csc_p_j, *pexsi_x_order);
        
        (*num_calls_pexsi)++;

    }else{

        update_x_pexsi(*a, *nnz, *dim, *pexsi_x_order, *pexsi_x);
        //print_x_pexsi(*pexsi_x, *nnz, *dim); // Dereference pexsi_x here
        
    }

    return 0; // Return 0 to indicate success
}

/*
    double* a = NULL;
    int* ia = NULL;
    int* p_ja = NULL;
    int nnz;
    int dim;

    int* sympack_csc_i = NULL;
    int* sympack_csc_p_j = NULL;
    int* sympack_x_order = NULL;
    int num_calls_sympack = 0;
    double* sympack_x = (double*)malloc(nnz * sizeof(double));


    int* pexsi_csc_i = NULL;
    int* pexsi_csc_p_j = NULL;
    int* pexsi_x_order = NULL;
    int num_calls_pexsi = 0;
    double* pexsi_x = (double*)malloc((2*nnz-dim) * sizeof(double));
    bool sympack_interface;
    bool pexsi_iterface;

*/


void wrapper_get_upper_csr_bin_c(
    const char* filename, 
    double** a, 
    int** ia, 
    int** p_ja, 
    int* nnz,
    int* dim, 
    bool zero_based,
    bool get_b, 
    bool get_x,
    double** b_rhs, double** sol_x, int* num_rhs) {

    FILE* file_bin = fopen(filename, "rb");
    if (file_bin == NULL) {
        perror("Error opening the file");
        return;
    }

    fread(dim, sizeof(int), 1, file_bin);
    fread(nnz, sizeof(int), 1, file_bin);
    if(get_b) fread(num_rhs, sizeof(int), 1, file_bin);

    *a = (double*)realloc(*a, *nnz * sizeof(double));
    *ia = (int*)realloc(*ia, *nnz * sizeof(int));
    *p_ja = (int*)realloc(*p_ja, (*dim + 1) * sizeof(int));
    if(get_b) *b_rhs = (double*)realloc(*b_rhs, (*num_rhs)*(*nnz) * sizeof(double));
    if(get_x) *sol_x = (double*)realloc(*sol_x, *nnz * sizeof(double));

    fread(*ia, sizeof(int), *nnz, file_bin);       // row indices
    fread(*p_ja, sizeof(int), *dim + 1, file_bin); // column pointers
    fread(*a, sizeof(double), *nnz, file_bin);     // values
    if(get_x) fread(*sol_x, sizeof(double), *nnz, file_bin);     // values
    if(get_b) fread(*b_rhs, sizeof(double), (*num_rhs)*(*nnz), file_bin);     // values
    fclose(file_bin);

    if (zero_based && (*ia)[0] == 1) {
        for (int i = 0; i < *nnz; i++) {
            (*ia)[i]--;
        }
        for (int i = 0; i < *dim + 1; i++) {
            (*p_ja)[i]--;
        }
    } else if (!zero_based && (*ia)[0] == 0) {
        for (int i = 0; i < *nnz; i++) {
            (*ia)[i]++;
        }
        for (int i = 0; i < *dim + 1; i++) {
            (*p_ja)[i]++;
        }
    }
}

void get_upper_csr_bin_c_args(int argc, char **argv, double** a, int** ia, int** p_ja, int* nnz, int* dim, bool get_b, bool get_x, double** b_rhs, double** sol_x, int* num_rhs) {

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input file name>" << std::endl;
        exit(1);  // Use exit to terminate the program or return an error code if you wish
    }
    wrapper_get_upper_csr_bin_c(argv[1], a, ia, p_ja, nnz, dim, false, get_b, get_x, b_rhs, sol_x, num_rhs);

}

void get_upper_csr_bin_c_filename(const char *filename, double** a, int** ia, int** p_ja, int* nnz, int* dim, bool get_b, bool get_x, double** b_rhs, double** sol_x, int* num_rhs) {

    wrapper_get_upper_csr_bin_c(filename, a, ia, p_ja, nnz, dim, false, get_b, get_x, b_rhs, sol_x, num_rhs);

}


