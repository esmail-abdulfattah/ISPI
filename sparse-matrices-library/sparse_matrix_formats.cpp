#include "sparse_matrix_formats.hpp"
#include <algorithm>  // for std::sort
#include <numeric>    // for std::iota
#include <cstdio>     // for FILE and associated functions
#include <iostream>   // for std::cout, std::endl
#include <iomanip>  // Needed for std::setprecision and std::fixed
#include <stdlib.h>
#include <ostream>
#include <vector>
#include <algorithm>
#include <cstring> // for std::memcpy
#include <fstream> // Include this header for file stream operations
#include <cmath> // Include this header for mathematical operations

struct NonZeroEntry {
    int row;
    int col;
    double value;
    size_t orig_idx;

    bool operator<(const NonZeroEntry &other) const {
        if (col != other.col) {
            return col < other.col;
        }
        return row < other.row;
    }
};

void SymmetricSparseMatrixFormats::readFromFile_upper_csr_bin(const char *filename) {

    FILE* file_bin = fopen(filename, "rb");
    if (file_bin == NULL) {
        perror("Error opening the file");
        return;
    }

    int nrows, nnz;
    fread(&nrows, sizeof(int), 1, file_bin);
    fread(&nnz, sizeof(int), 1, file_bin);
    
    csr_x.resize(nnz);
    csr_i.resize(nnz);
    csr_p_j.resize(nrows + 1);

    upper_nnz = nnz - nrows;
    dim = nrows;

    fread(csr_i.data(), sizeof(int), nnz, file_bin);       // row indices
    fread(csr_p_j.data(), sizeof(int), nrows + 1, file_bin); // column pointers
    fread(csr_x.data(), sizeof(double), nnz, file_bin);    // values
    tick_readFromFile_upper_csr_bin = true;
    fclose(file_bin);


    if (zero_based && csr_i[0] == 1) {
        for (auto& val : csr_i) {
            val--;
        }
        for (auto& val : csr_p_j) {
            val--;
        }
    }else if(!zero_based && csr_i[0] == 0){
        for (auto& val : csr_i) {
            val++;
        }
        for (auto& val : csr_p_j) {
            val++;
        }
    }
}

void SymmetricSparseMatrixFormats::print_upper_csr() {

    if(tick_readFromFile_upper_csr_bin){
        std::cout << "csr_i (Row indices): ";
        for (const auto &index : csr_i) {
            std::cout << index << " ";
        }
        std::cout << std::endl;

        std::cout << "csr_p_j (Column pointers): ";
        for (const auto &pointer : csr_p_j) {
            std::cout << pointer << " ";
        }
        std::cout << std::endl;

        std::cout << "csr_x (Values): ";
        for (const auto &value : csr_x) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }else{
        std::cerr << "Error: The matrix hasn't been successfully read from the file. Please call readFromFile_upper_csr_bin() first." << std::endl;
    }
    
}

void SymmetricSparseMatrixFormats::get_full_csr_column_indices_from_csr_column_pointers() {
    
    int n_cols = csr_p_j.size() - 1;
    
    int start = 1;
    if(zero_based) start = 0;
    for (int index_j = start; index_j <= n_cols; ++index_j) {  // Start loop from 1
        for (int count = 0; count < (csr_p_j[index_j] - csr_p_j[index_j - 1]); ++count) {
            //std::cout << index_j << std::endl;
            csr_j.push_back(index_j);
        }
    }
}

// Define a function that converts the upper CSR representation of a symmetric matrix
// to the lower CSC representation.
int SymmetricSparseMatrixFormats::get_lower_csc_from_upper_csr(bool getx) {
    
    // Check if the conversion has been performed previously to avoid redundancy.

    if(false){
    if(!tick_get_lower_csc_from_upper_csr){

        // Ensure the CSR column indices are complete.
        if (csr_j.empty()) get_full_csr_column_indices_from_csr_column_pointers();

        // Verify that the vectors for row indices, column indices, and values have the same size.
        if (csr_i.size() != csr_j.size() || csr_i.size() != csr_x.size()) {
            std::cerr << "The vectors have inconsistent sizes!" << std::endl;
            return 1; // Return error code for inconsistent sizes.
        }

        // Create a temporary vector to hold non-zero entries.
        std::vector<NonZeroEntry> entries;

        // Populate the entries with data from the CSR format.
        for (size_t k = 0; k < csr_i.size(); ++k) {
            entries.push_back({csr_i[k], csr_j[k], csr_x[k], k});
        }

        // Sort the entries to prepare for CSC format.
        std::sort(entries.begin(), entries.end());
    
        // Initialize the first column index.
        int currentCol = entries.empty() ? -1 : entries[0].col; 
        // The CSC 'p' vector starts with a 1 due to 1-based indexing.
        csc_p_j.push_back(1);  

        // Convert the sorted entries to CSC format.
        for (const auto &entry : entries) {
            csc_i.push_back(entry.row); // Store the row indices.
            if (getx) csc_x.push_back(entry.value); // Conditionally store the values.
            csc_sorted_order_.push_back(entry.orig_idx); // Store the original indices of the entries.

            // When a new column is encountered, add the current size of `csc_i` to `csc_p_j`.
            if (entry.col != currentCol) {
                csc_p_j.push_back(csc_i.size());
                currentCol = entry.col;
            }
        }
        // Finalize the `csc_p_j` vector with the size of `csc_i` plus one.
        csc_p_j.push_back(csc_i.size() + 1); 

    }
    // Mark the conversion as done.
    tick_get_lower_csc_from_upper_csr = true;
    // If everything is successful, return 0.
    }

    if(!tick_get_lower_csc_from_upper_csr){

        if (csr_j.empty()) get_full_csr_column_indices_from_csr_column_pointers();

        // Reserve maximum possible space
        csc_i = csr_j;
        csc_j = csr_i;
        if(getx) csc_x = csr_x;

        // Sort by columns first and then by rows
        csc_sorted_order_.resize(csc_i.size());
        std::iota(csc_sorted_order_.begin(), csc_sorted_order_.end(), 0); // Fill it with 0 to dim-1.
        std::sort(csc_sorted_order_.begin(), csc_sorted_order_.end(),
        [this](size_t a, size_t b) {
            return csc_j[a] < csc_j[b] || (csc_j[a] == csc_j[b] && csc_i[a] < csc_i[b]);
        });


        std::vector<int> new_csc_i(csc_i.size());
        std::vector<int> new_csc_j(csc_j.size());
        std::vector<double> new_csc_x(csc_x.size());


        for (size_t i = 0; i < csc_sorted_order_.size(); ++i) {
            new_csc_i[i] = csc_i[csc_sorted_order_[i]];
            new_csc_j[i] = csc_j[csc_sorted_order_[i]];
            if(getx) new_csc_x[i] = csc_x[csc_sorted_order_[i]];
        }

        csc_i = new_csc_i;
        csc_j = new_csc_j;
        if(getx) csc_x = new_csc_x; 

        csc_p_j.clear();
        csc_p_j.resize(csc_j.back() + 1, 0);  // Adjusted resize length
        for (const auto& j_val : csc_j) {
            csc_p_j[j_val]++;  // Adjusted indexing
        }
        for (size_t idx = 1; idx < csc_p_j.size(); ++idx) {
            csc_p_j[idx] += csc_p_j[idx - 1];
        }

        if (!zero_based) {
        for (auto& val : csc_p_j) {
            val += 1;
        }
    }
    
    tick_get_lower_csc_from_upper_csr = true;
    
    }
    return 0;
}

void SymmetricSparseMatrixFormats::print_lower_csc() {

    if(tick_get_lower_csc_from_upper_csr){
        std::cout << "csc_i: ";
        for (const auto &val : csc_i) std::cout << val << " ";
        std::cout << std::endl;

        std::cout << "csc_p_j: ";
        for (const auto &val : csc_p_j) std::cout << val << " ";
        std::cout << std::endl;

        std::cout << "csc_x: ";
        for (const auto &val : csc_x) std::cout << val << " ";
        std::cout << std::endl;

        std::vector<double> reordered_x(csc_x.size());
        for (size_t i = 0; i < csc_x.size(); ++i) {
            reordered_x[i] = csr_x[csc_sorted_order_[i]];
        }

        std::cout << "reordered_x: ";
        for (const auto &val : reordered_x) std::cout << val << " ";
        std::cout << std::endl;

        std::cout << "Order from x_lower to csc_x: ";
        for (const auto &idx : csc_sorted_order_) std::cout << idx << " ";
        std::cout << std::endl;

    }else{
        std::cerr << "Error: The matrix hasn't been successfully read from the file. Please call readFromFile_upper_csr_bin() first." << std::endl;
    }
    
}


// Function to assign CSC vectors to global pointers
void SymmetricSparseMatrixFormats::exportToGlobalLowerCSCPointers(int*& global_csc_i, 
                                int*& global_csc_p_j, 
                                int*& global_csc_x_order) {
    // Allocate memory for the global pointers
    global_csc_i = new int[csc_i.size()];
    global_csc_p_j = new int[csc_p_j.size()];
    global_csc_x_order = new int[csc_sorted_order_.size()];

    // Copy the data from the class vectors to the allocated arrays
    std::copy(csc_i.begin(), csc_i.end(), global_csc_i);
    std::copy(csc_p_j.begin(), csc_p_j.end(), global_csc_p_j);
    std::copy(csc_sorted_order_.begin(), csc_sorted_order_.end(), global_csc_x_order);
}

// Function to assign CSC vectors to global pointers
void SymmetricSparseMatrixFormats::exportToGlobalFullCSCPointers(int*& global_csc_i, 
                                int*& global_csc_p_j, 
                                int*& global_csc_x_order) {
    // Allocate memory for the global pointers
    global_csc_i = new int[full_i.size()];
    global_csc_p_j = new int[full_p_j.size()];
    global_csc_x_order = new int[full_sorted_order_.size()];

    // Copy the data from the class vectors to the allocated arrays
    std::copy(full_i.begin(), full_i.end(), global_csc_i);
    std::copy(full_p_j.begin(), full_p_j.end(), global_csc_p_j);
    std::copy(full_sorted_order_.begin(), full_sorted_order_.end(), global_csc_x_order);
}

void SymmetricSparseMatrixFormats::initialize_from_upper_triangular(bool getx) {
    
    if(false){
        if (csr_j.empty()) get_full_csr_column_indices_from_csr_column_pointers();

        full_i = csr_i;
        full_j = csr_j;
        full_x = csr_x;  // This will change each time, but initializing once
            
        // Reserve maximum possible space
        full_i.reserve(2 * upper_nnz + dim);
        full_j.reserve(2 * upper_nnz + dim);
        full_x.reserve(2 * upper_nnz + dim);

        std::vector<int> indices;
        indices.reserve(upper_nnz + dim);

        // Create mirrored elements and append directly to full_i, full_j, and full_x
        for (size_t idx = 0; idx < csr_i.size(); ++idx) {
            if (csr_i[idx] != csr_j[idx]) {
                full_i.push_back(csr_j[idx]);
                full_j.push_back(csr_i[idx]);
                full_x.push_back(csr_x[idx]);  // This line could be skipped if csr_x is going to be replaced every time
                indices.push_back(idx);
            }
        }

        full_sorted_order_.resize(csr_i.size());
        std::iota(full_sorted_order_.begin(), full_sorted_order_.end(), 0); // Fill it with 0 to dim-1.
        full_sorted_order_.insert(full_sorted_order_.end(), indices.begin(), indices.end());

        std::sort(full_sorted_order_.begin(), full_sorted_order_.end(),
                [this](size_t a, size_t b) {
                    return full_j[a] < full_j[b] || (full_j[a] == full_j[b] && full_i[a] < full_i[b]);
                });

    }else{
        if (csr_j.empty()) get_full_csr_column_indices_from_csr_column_pointers();

        full_i = csr_i;
        full_j = csr_j;
        if(getx) full_x = csr_x;  

        // Reserve maximum possible space
        full_i.reserve(2 * upper_nnz + dim);
        full_j.reserve(2 * upper_nnz + dim);
        if(getx) full_x.reserve(2 * upper_nnz + dim);
        std::vector<int> indices;
        indices.reserve(upper_nnz);

        // Create mirrored elements and append directly to full_i, full_j, and full_x
        for (size_t idx = 0; idx < csr_i.size(); ++idx) {
            if (csr_i[idx] != csr_j[idx]) {
                full_i.push_back(csr_j[idx]);
                full_j.push_back(csr_i[idx]);
                if(getx) full_x.push_back(csr_x[idx]);  // This line could be skipped if csr_x is going to be replaced every time
                indices.push_back(idx);
            }
        }

        // Sort by columns first and then by rows
        full_sorted_order_.resize(csr_i.size());
        std::iota(full_sorted_order_.begin(), full_sorted_order_.end(), 0); // Fill it with 0 to dim-1.
        full_sorted_order_.insert(full_sorted_order_.end(), indices.begin(), indices.end());


        std::vector<size_t> sorted_indices(full_i.size());
        std::iota(sorted_indices.begin(), sorted_indices.end(), 0);

        // Sort the indices based on the sorting criteria used for full_sorted_order_
        std::sort(sorted_indices.begin(), sorted_indices.end(),
        [this](size_t a, size_t b) {
            return full_j[a] < full_j[b] || (full_j[a] == full_j[b] && full_i[a] < full_i[b]);
        });

        std::vector<int> new_order(full_sorted_order_.size());
        for (int i = 0; i < sorted_indices.size(); ++i) {
            new_order[i] = full_sorted_order_[sorted_indices[i]];
        }

        // Replace full_sorted_order_ with the new order
        full_sorted_order_ = new_order;

        // If you want to reorder other vectors like full_i, full_j, and full_x, you can do it similarly:
        std::vector<int> new_full_i(full_i.size());
        std::vector<int> new_full_j(full_j.size());
        std::vector<double> new_full_x(full_x.size());

        for (int i = 0; i < sorted_indices.size(); ++i) {
            new_full_i[i] = full_i[sorted_indices[i]];
            new_full_j[i] = full_j[sorted_indices[i]];
            if(getx) new_full_x[i] = csr_x[full_sorted_order_[i]];
        }

        // Now replace the old vectors with the new ordered vectors
        full_i = new_full_i;
        full_j = new_full_j;
        if(getx) full_x = new_full_x;

        /*
        std::cout << std::endl;
        std::cout << "full_i: ";
        for (const auto &pointer : full_i) {
            std::cout << pointer << " ";
        }
        std::cout << std::endl;


        std::cout << std::endl;
        std::cout << "full_j: ";
        for (const auto &pointer : full_j) {
            std::cout << pointer << " ";
        }
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << "full_x: ";
        for (const auto &pointer : full_x) {
            std::cout << pointer << " ";
        }
        std::cout << std::endl;


        exit(0);
        */

      }
      

}

void SymmetricSparseMatrixFormats::update_x(std::vector<double>& values) {

    csr_x = values;
    
    full_x = csr_x;
    for (size_t idx = 0; idx < csr_i.size(); ++idx) {
        if (csr_i[idx] != csr_j[idx]) {
            full_x.push_back(csr_x[idx]);  // This line could be skipped if csr_x is going to be replaced every time
        }
    }
/*
    if(!values.empty()){

        //std::cout << "values size: " << values.size() << std::endl;
        //std::cout << "full_sorted_order_ size: " << full_sorted_order_.size() << std::endl;

        for (size_t idx = 0; idx < full_sorted_order_.size(); ++idx) {
            full_x[idx] = values[full_sorted_order_[idx]];
        }

        //for (const auto &index : full_sorted_order_) {
        //    std::cout << index << " ";
        //}
        //std::cout << std::endl;
    }
*/

}

void SymmetricSparseMatrixFormats::reconstruct_column_pointers(){

    // Re-construct the column pointer vector
    full_p_j.clear();
    full_p_j.resize(full_j.back() + 1, 0);  // Adjusted resize length
    for (const auto& j_val : full_j) {
        full_p_j[j_val]++;  // Adjusted indexing
    }
    for (size_t idx = 1; idx < full_p_j.size(); ++idx) {
        full_p_j[idx] += full_p_j[idx - 1];
    }

    if (!zero_based) {
        for (auto& val : full_p_j) {
            val += 1;
        }
    }

}

void SymmetricSparseMatrixFormats::get_full_csr_from_upper_csr(bool getx){

    if(!tick_get_full_csc_from_upper_csr){
        initialize_from_upper_triangular(getx);
        if(getx) update_x(csr_x);
        reconstruct_column_pointers();
    }
    tick_get_full_csc_from_upper_csr = true;
}

void SymmetricSparseMatrixFormats::get_full_csc_from_upper_csr(bool getx){

    get_full_csr_from_upper_csr(getx);

}

void SymmetricSparseMatrixFormats::print_full_csc() {

    if(tick_get_full_csc_from_upper_csr){
        std::cout << "full_i (Row indices): ";
        for (const auto &index : full_i) {
            std::cout << index << " ";
        }
        std::cout << std::endl;

        std::cout << "full_p_j (Column Pointers): ";
        for (const auto &pointer : full_p_j) {
            std::cout << pointer << " ";
        }
        std::cout << std::endl;

        std::cout << "full_x (Values): ";
        for (const auto &value : full_x) {
            std::cout << value << " ";
        }
        std::cout << std::endl;

        std::cout << std::endl;
        std::cout << "full_sorted_order_: ";
        for (const auto &value : full_sorted_order_) {
            std::cout << value << " ";
        }
        std::cout << std::endl;

    }else{
        std::cerr << "Error: The matrix hasn't been successfully read from the file. Please call readFromFile_upper_csr_bin() first." << std::endl;
    }
    
}

void SymmetricSparseMatrixFormats::print_full_csr() {

    print_full_csc();
}

int SymmetricSparseMatrixFormats::export_rb_format(const std::string &inputFileName) {

    get_lower_csc_from_upper_csr();

    // Definition of argv as an empty array.
    const char* argv[] = {};
    int nnz = upper_nnz + dim;

    std::string outputFileName = inputFileName.substr(0, inputFileName.find_last_of('.')) + ".rb";

    const char* file = outputFileName.c_str();
    FILE* fid = fopen(file, "wt+");
    if (!fid) {
        perror("Can't open file for writing");
        return -1;
    }

    const char* key = sizeof(argv) > 4 ? argv[4] : "key";
    const char* title = sizeof(argv) > 3 ? argv[3] : "title";
    const char* type = sizeof(argv) > 5 ? argv[5] : "RUA";
    int ifmt = sizeof(argv) > 6 ? atoi(argv[6]) : 14;
    int job = sizeof(argv) > 7 ? atoi(argv[7]) : 2;

    // Rest of your code, modified as per your instructions ...

    int ncol = dim; // Assuming square matrix
    int len = ceil(log10(0.1 + upper_nnz + dim + 1)) + 1;
    int nperline = std::min((int)floor(80.0/len), ncol+1);
    int ptr_len = len;
    int ptr_nperline = nperline;
    int ptrcrd = floor(ncol / nperline) + 1;
    char ptrfmt[100];
    sprintf(ptrfmt, "(%dI%d)", nperline, len);
    
    // Compute row index format
    nperline = std::min((int)floor(80.0/len), upper_nnz + dim);
    int ind_len = len;
    int ind_nperline = nperline;
    int indcrd = floor((nnz-1) / nperline) + 1;
    char indfmt[100];
    sprintf(indfmt, "(%dI%d)", nperline, len);

    // The rest of your code
    // Compute values and rhs format (same format for both)
    int valcrd = 0;
    int rhscrd = 0;
    std::string c_valfmt;
    std::string valfmt;

    if (job > 1) {
        int ihead = 0;
        if (ifmt >= 100) {
            ihead = ifmt / 100;
            ifmt = ifmt - 100 * ihead;
        }

        len = ifmt + (ifmt >= 100 ? ihead + 2 : 7);
        nperline = floor(80.0 / len);
        int c_len = len;

        c_valfmt = "%" + std::to_string(c_len) + "." + std::to_string(ifmt) + (ifmt >= 100 ? "f" : "E");
        valfmt = std::to_string(nperline) + (ifmt >= 100 ? "F" : "E") + std::to_string(len) + "." + std::to_string(ifmt);
        valcrd = floor((nnz - 1) / nperline) + 1;
        valfmt = "(" + valfmt + ")";
    }


    int nrhs = job - 2;
    if (nrhs >= 1) {
        rhscrd = floor((nrhs * dim - 1) / nperline) + 1;
    }

    int totcrd = ptrcrd + indcrd + valcrd + rhscrd;

    // Line 1
    std::string t = title;
    int m = t.size();
    for (int i = m + 1; i <= 72; ++i) {
        t += " ";
    }
    fprintf(fid, "%72s", t.c_str());
    t = key;
    m = t.size();
    for (int i = m + 1; i <= 8; ++i) {
        t += " ";
    }
    fprintf(fid, "%8s\n", t.c_str());

    // Line 2
    fprintf(fid, "%14d%14d%14d%14d%14d\n", totcrd, ptrcrd, indcrd, valcrd, rhscrd);

    // Line 3
    t = type;
    m = t.size();
    for (int i = m + 1; i <= 14; ++i) {
        t += " ";
    }
    fprintf(fid, "%14s", t.c_str());
    fprintf(fid, "%14i%14i%14i%14i\n", dim, ncol, nnz, nrhs); // Assuming nnzeros is nnz and nrow is nrows

    // Line 4
    t = ptrfmt; 
    m = t.size();
    for (int i = m + 1; i <= 16; ++i) {
        t += " ";
    }
    fprintf(fid, "%16s", t.c_str());

    t = indfmt;
    m = t.size();
    for (int i = m + 1; i <= 16; ++i) {
        t += " ";
    }
    fprintf(fid, "%16s", t.c_str());

    t = valfmt;
    m = t.size();
    for (int i = m + 1; i <= 20; ++i) {
        t += " ";
    }
    fprintf(fid, "%20s", t.c_str());

    t = valfmt;  // Updated this line to reuse valfmt for the RHS vector as well
    m = t.size();
    for (int i = m + 1; i <= 20; ++i) {
        t += " ";
    }
    fprintf(fid, "%20s\n", t.c_str());


    // column pointers
    t = " ";
    for (int j = 1; j <= ptr_nperline; ++j) {
        t += "%" + std::to_string(ptr_len) + "d ";
    }
    t += "\n";
    for (int i = 0; i < ncol + 1; i += ptr_nperline) {
        int limit = std::min(i + ptr_nperline, ncol + 1);
        for (int j = i; j < limit; ++j) {
            std::string format_string = "%" + std::to_string(ptr_len) + "d";
            fprintf(fid, format_string.c_str(), csc_p_j[j]); // using csr_p_j
        }
        fprintf(fid, "\n");
    }

    std::string formatStr = "%" + std::to_string(ind_len) + "d";
    for (int i = 0; i < nnz; ++i) {
        if (i > 0 && i % ind_nperline == 0) {
            fprintf(fid, "\n");
        }
        fprintf(fid, formatStr.c_str(), csc_i[i]); // using i
    }

    fprintf(fid, "\n");
    if (job >= 2) {
        for (int i = 0; i < nnz; ++i) {
            fprintf(fid, c_valfmt.c_str(), csc_x[i]); // using csr_x
            if ((i + 1) % nperline == 0) {
                fprintf(fid, "\n");
            }
        }
        if (nnz % nperline != 0) {
            fprintf(fid, "\n");
        }
    }
    
    fclose(fid);

    return 0;
}



/*
void SymmetricSparseMatrixFormats::csc_lower_triangular_column_format(){

    std::vector<int> lower_tri_i;
    std::vector<int> lower_tri_j;
    std::vector<double> lower_tri_x;

    int nnz = dim + upper_nnz;
    p_i_csc.resize(dim + 1, 0);
    j_csc.resize(nnz, 0);
    x_csc.resize(nnz, 0);

    // Filter out upper triangular elements
    for(size_t index = 0; index < full_i.size(); ++index) {
        if(full_j[index] <= full_i[index]) {
            lower_tri_i.push_back(full_i[index]);
            lower_tri_j.push_back(full_j[index]);
            lower_tri_x.push_back(full_x[index]);
        }
    }

    // Populate ja and values using lower triangular data
    j_csc = lower_tri_i;
    x_csc = lower_tri_x;

    // Populate ia using column indices from lower_tri_j
    for(const int& col : lower_tri_j) {
        p_i_csc[col + 1]++;
    }

    for(int i = 1; i < dim; ++i) {
        p_i_csc[i + 1] += p_i_csc[i];
    }
}

int SymmetricSparseMatrixFormats::sort_columns() {

    int n = i.size(); // Assuming 'i' is the replacement for 'rows'

    // Step 1: Create a vector of indices [0, 1, 2, ...]
    std::vector<int> indices(n);
    for (int idx = 0; idx < n; ++idx) {
        indices[idx] = idx;
    }
    
    // Step 2: Sort the indices vector based on the values in the i vector

    std::sort(indices.begin(), indices.end(), [this](int a, int b) {
        return i[a] < i[b]; 
    });


    // Step 3: Create new vectors to hold the sorted rows and columns
    std::vector<int> sortedI(n);
    std::vector<int> sortedP_j(n);
    std::vector<double> sortedX(n);

    // Use the sorted indices to rearrange i, csr_p_j, and csr_x
    for (int idx = 0; idx < n; ++idx) {
        sortedI[idx] = i[indices[idx]];
        sortedP_j[idx] = csr_p_j[indices[idx]];
        sortedX[idx] = csr_x[indices[idx]];
    }

    // Now, move the sorted data back to the original vectors
    i = sortedI;
    csr_p_j = sortedP_j;
    csr_x = sortedX;

    return 0;
}


*/