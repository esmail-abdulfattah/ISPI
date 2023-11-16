/**
 * @file sparse_matrix_formats.hpp
 * @author Esmail Abdul Fattah
 * @date Date Created Octobor 2023
 * 
 * Description: This header defines the SymmetricSparseMatrixFormats class and associated
 * methods for handling sparse matrix data in different format.
 * 
 * License (optional): This is free and unencumbered software released into the public domain.
 * Anyone is free to copy, modify, publish, use, compile, sell, or
 * distribute this software, either in source code form or as a compiled
 * binary, non-commercial way.
 * 
 */



#ifndef SPARSE_MATRIX_FORMATS_H
#define SPARSE_MATRIX_FORMATS_H

#include <vector>
#include <string>


class SymmetricSparseMatrixFormats {
private:

    int upper_nnz = 0, dim = 0;

    //upper csr format:
    std::vector<int> csr_i;    // Row indices
    std::vector<int> csr_j;    // Column indices
    std::vector<double> csr_x; // Values
    std::vector<int> csr_p_j;    // Column pointers
    bool tick_readFromFile_upper_csr_bin = false;

    //get lower csc from upper csr
    std::vector<int> csc_i;
    std::vector<int> csc_j;
    std::vector<int> csc_p_j; 
    std::vector<double> csc_x;
    std::vector<int> csc_sorted_order_; // Keeps track of the original order
    bool tick_get_lower_csc_from_upper_csr = false;

    //full csc/csr format
    std::vector<int> full_i;    
    std::vector<int> full_j;
    std::vector<double> full_x; 
    std::vector<int> full_p_j;    
    std::vector<int> full_sorted_order_;
    bool tick_get_full_csc_from_upper_csr = false;

public:

    bool zero_based = false;

    // Constructor
    SymmetricSparseMatrixFormats() = default; // Add this line
    SymmetricSparseMatrixFormats(double* a, int* ia, int* p_ja, int nnz_val, int dim_val) : upper_nnz(nnz_val - dim_val), dim(dim_val) {
        int nnz = nnz_val;
        csr_i.assign(ia, ia + nnz);
        csr_x.assign(a, a + nnz);
        csr_p_j.assign(p_ja, p_ja + dim + 1);
        tick_readFromFile_upper_csr_bin = true;
    }

    // Constructor using vectors
    SymmetricSparseMatrixFormats(const std::vector<double>& a, 
                          const std::vector<int>& ia, 
                          const std::vector<int>& p_ja, 
                          int upper_nnz_val, 
                          int dim_val) 
        : csr_x(a), csr_i(ia), csr_p_j(p_ja), upper_nnz(upper_nnz_val), dim(dim_val) 
    {
        tick_readFromFile_upper_csr_bin = true;
    }

    //csr format:
    void readFromFile_upper_csr_bin(const char *filename);
    void print_upper_csr();

    //rb format:
    int export_rb_format(const std::string &inputFileName);
    //int sort_columns();
    //void csc_lower_triangular_column_format();

    //convert to full matrix
    void get_full_csr_column_indices_from_csr_column_pointers();
    void initialize_from_upper_triangular(bool getx = false);
    void update_x(std::vector<double>& values);
    void reconstruct_column_pointers();
    void get_full_csr_from_upper_csr(bool getx = false);
    void get_full_csc_from_upper_csr(bool getx = false);
    void print_full_csc();
    void print_full_csr();

    int get_lower_csc_from_upper_csr(bool getx = false);
    void print_lower_csc();

    //export functions:
    //Function to assign CSC vectors to global pointers
    void exportToGlobalLowerCSCPointers(int*& global_csc_i, int*& global_csc_p_j, int*& global_csc_x_order);                               
    void exportToGlobalFullCSCPointers(int*& global_csc_i, int*& global_csc_p_j, int*& global_csc_x_order);


};


#endif // SPARSE_MATRIX_FORMATS_H
