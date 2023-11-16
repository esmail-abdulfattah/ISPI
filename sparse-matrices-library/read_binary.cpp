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

int sort_columns(int *rows, int *columns, int n, int *sortedColumns, double *values);
void csc_lower_triangular_column_format(int nrows, int nnz, int *rows, int *columns, int *ia, int *ja);
void csr_lower_triangular_row_format(int nrows, int nnz, int *rows, int *columns, int *ia, int *ja, double *values);


int main_1() {

    FILE *file = fopen("Q.bin", "rb");

    if (file == NULL) {
        perror("Error opening the file");
        return 1;
    }

    int nrows, nnz;

    fread(&nrows, sizeof(int), 1, file);
    fread(&nnz, sizeof(int), 1, file);

    printf("%d %d\n", nrows, nnz);
    
    int *rows = (int *)malloc(nnz * sizeof(int));
    int *columns = (int *)malloc(nnz * sizeof(int));
    double *values = (double *)malloc(nnz * sizeof(double));

    fread(rows, sizeof(int), nnz, file); // row indices
    fread(columns, sizeof(int), nnz, file); // column indices
    fread(values, sizeof(double), nnz, file); // values

    //std::cout << "Rows and columns" << std::endl;
    //for(int i=0; i < nnz; i++) std::cout << rows[i] << ", "; std::cout << std::endl; //csc formatted
    //for(int i=0; i < nnz; i++) std::cout << columns[i] << ", "; std::cout << std::endl; //row index!
    
    //--------------------------------------------------------------------------------------------------------------------------------------------------
    //csc format:
    //--------------------------------------------------------------------------------------------------------------------------------------------------
    int *ia_csc = (int *)malloc(nnz * sizeof(int));
    int *ja_csc = (int *)malloc((nrows + 1) * sizeof(int));
    csc_lower_triangular_column_format(nrows, nnz, rows, columns, ia_csc, ja_csc);

    //std::cout << "csc format" << std::endl;
    //for(int i=0; i <=nrows; i++) std::cout << ja_csc[i] << ", "; std::cout << std::endl; 
    //for(int i=0; i < nnz; i++) std::cout << ia_csc[i] << ", "; std::cout << std::endl; 
    //for(int i=0; i < nnz; i++) std::cout << values[i] << ", "; std::cout << std::endl; 
    //--------------------------------------------------------------------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------------------------------------------------------------------
    //csr format:
    //--------------------------------------------------------------------------------------------------------------------------------------------------
    int *ia_csr = (int *)malloc((nrows + 1) * sizeof(int));
    int *ja_csr = (int *)malloc(nnz * sizeof(int));

    csr_lower_triangular_row_format(nrows, nnz, rows, columns, ia_csr, ja_csr, values);

    //std::cout << "csr format" << std::endl;
    //for(int i=0; i <=nrows; i++) std::cout << ia_csr[i] << ", "; std::cout << std::endl; 
    //for(int i=0; i < nnz; i++) std::cout << ja_csr[i] << ", "; std::cout << std::endl; 
    //for(int i=0; i < nnz; i++) std::cout << values[i] << ", "; std::cout << std::endl;
    //--------------------------------------------------------------------------------------------------------------------------------------------------

    return 0;

}


int main(int argc, char* argv[]) {
    
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_filename> [other_arguments...]\n";
        return -1;
    }

    std::string inputFileName = argv[1]; // The input file name is the first argument
    std::string outputFileName = inputFileName.substr(0, inputFileName.find_last_of('.')) + ".rb";

    std::ofstream outFile(outputFileName);
    if (!outFile.is_open()) {
        std::perror("Can't open file for writing");
        return -1;
    }


    const char* file = outputFileName.c_str();
    FILE* fid = fopen(file, "wt+");
    if (!fid) {
        perror("Can't open file for writing");
        return -1;
    }

    const char* key = argc > 4 ? argv[4] : "key";
    const char* title = argc > 3 ? argv[3] : "title";
    const char* type = argc > 5 ? argv[5] : "RUA";
    int ifmt = argc > 6 ? atoi(argv[6]) : 14;
    int job = argc > 7 ? atoi(argv[7]) : 2;

    FILE* file_bin = fopen(inputFileName.c_str(), "rb");

    if (file_bin == NULL) {
        perror("Error opening the file");
        return 1;
    }

    int nrows, nnz;
    fread(&nrows, sizeof(int), 1, file_bin);
    fread(&nnz, sizeof(int), 1, file_bin);

    std::cout << nrows << " ------- " << nnz << std::endl;
    printf("%d %d\n", nrows, nnz);
    
    if(false){
        int* rows = (int*)malloc(nnz * sizeof(int));
        int* columns = (int*)malloc(nnz * sizeof(int));
        double* a = (double*)malloc(nnz * sizeof(double));

        fread(rows, sizeof(int), nnz, file_bin); // row indices
        fread(columns, sizeof(int), nnz, file_bin); // column indices
        fread(a, sizeof(double), nnz, file_bin); // values

        int* ia = (int*)malloc(nnz * sizeof(int));
        int* ja = (int*)malloc((nrows + 1) * sizeof(int));
        csc_lower_triangular_column_format(nrows, nnz, rows, columns, ia, ja);

        for(int i =0 ;i < (nrows + 1); i++) std::cout << ja[i] << ", "; std::cout << std::endl; std::cout << std::endl;
        for(int i =0 ;i < nnz; i++) std::cout << ia[i] << ", "; std::cout << std::endl; std::cout << std::endl;

    }

    double* a = (double*)malloc(nnz * sizeof(double));
    int* ia = (int*)malloc(nnz * sizeof(int));
    int* ja = (int*)malloc((nrows + 1) * sizeof(int));
    fread(ia, sizeof(int), nnz, file_bin); // row indices
    fread(ja, sizeof(int), (nrows + 1), file_bin); // column indices
    fread(a, sizeof(double), nnz, file_bin); // values

    //for(int i =0 ;i < nnz; i++) std::cout << ia[i] << ", "; std::cout << std::endl; std::cout << std::endl;
    //for(int i =0 ;i < (nrows + 1); i++) std::cout << ja[i] << ", "; std::cout << std::endl; std::cout << std::endl;

    int ncol = nrows; // Assuming square matrix
    int len = ceil(log10(0.1 + nnz + 1)) + 1;
    int nperline = std::min((int)floor(80.0/len), ncol+1);
    int ptr_len = len;
    int ptr_nperline = nperline;
    int ptrcrd = floor(ncol / nperline) + 1;
    char ptrfmt[100];
    sprintf(ptrfmt, "(%dI%d)", nperline, len);
    
    // Compute row index format
    nperline = std::min((int)floor(80.0/len), nnz);
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
        int nrow = nrows; // Assuming nrow to be nrows, adjust as necessary
        rhscrd = floor((nrhs * nrow - 1) / nperline) + 1;
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
    fprintf(fid, "%14i%14i%14i%14i\n", nrows, ncol, nnz, nrhs); // Assuming nnzeros is nnz and nrow is nrows

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
            fprintf(fid, format_string.c_str(), ja[j]);
        }
        fprintf(fid, "\n");
    }

    // Preparing the format string
    std::string formatStr = "%" + std::to_string(ind_len) + "d";

    // Loop to print all elements with ind_nperline elements per line
    for (int i = 0; i < nnz; ++i) {
        if (i > 0 && i % ind_nperline == 0) {
            fprintf(fid, "\n");
        }
        fprintf(fid, formatStr.c_str(), ia[i]);
    }

    // Add a newline if the last line isn't full
    if (nnz % nperline != 0) {
        fprintf(fid, "\n");
    }

    // Check if a newline is needed at the end
    if (nnz==ind_nperline) {
        fprintf(fid, "\n");
    }

    // Checking the initialization of the 'a' array (Add this before your code to print the values)
    //for (int i = 0; i < nnz; ++i) {
    //    std::cout << "a[" << i << "] = " << a[i] << std::endl;
    //}

    // numerical values of nonzero elements of the matrix
    if (job >= 2) {
        for (int i = 0; i < nnz; ++i) {
            fprintf(fid, c_valfmt.c_str(), a[i]);
            if ((i + 1) % nperline == 0) {
                fprintf(fid, "\n");
            }
        }

        // Add a newline if the last line isn't full
        if (nnz % nperline != 0) {
            fprintf(fid, "\n");
        }
    }
    
    fclose(fid);
    fclose(file_bin);
    free(a);
    free(ia);
    free(ja);

    return 0;
}


int main_2(int argc, char* argv[]) {

    const char* file = "Q.rb";
    FILE* fid = fopen(file, "wt+");
    if (!fid) {
        perror("Can't open file for writing");
        return -1;
    }

    const char* key = argc > 4 ? argv[4] : "key";
    const char* title = argc > 3 ? argv[3] : "title";
    const char* type = argc > 5 ? argv[5] : "RUA";
    int ifmt = argc > 6 ? atoi(argv[6]) : 14;
    int job = argc > 7 ? atoi(argv[7]) : 2;

    FILE* file_bin = fopen("Q.bin", "rb");

    if (file_bin == NULL) {
        perror("Error opening the file");
        return 1;
    }

    int nrows, nnz;
    fread(&nrows, sizeof(int), 1, file_bin);
    fread(&nnz, sizeof(int), 1, file_bin);

    std::cout << nrows << " ------- " << nnz << std::endl;
    printf("%d %d\n", nrows, nnz);
    
    if(false){
        int* rows = (int*)malloc(nnz * sizeof(int));
        int* columns = (int*)malloc(nnz * sizeof(int));
        double* a = (double*)malloc(nnz * sizeof(double));

        fread(rows, sizeof(int), nnz, file_bin); // row indices
        fread(columns, sizeof(int), nnz, file_bin); // column indices
        fread(a, sizeof(double), nnz, file_bin); // values

        int* ia = (int*)malloc(nnz * sizeof(int));
        int* ja = (int*)malloc((nrows + 1) * sizeof(int));
        csc_lower_triangular_column_format(nrows, nnz, rows, columns, ia, ja);

        for(int i =0 ;i < (nrows + 1); i++) std::cout << ja[i] << ", "; std::cout << std::endl; std::cout << std::endl;
        for(int i =0 ;i < nnz; i++) std::cout << ia[i] << ", "; std::cout << std::endl; std::cout << std::endl;

    }

    double* a = (double*)malloc(nnz * sizeof(double));
    int* ia = (int*)malloc(nnz * sizeof(int));
    int* ja = (int*)malloc((nrows + 1) * sizeof(int));
    fread(ia, sizeof(int), nnz, file_bin); // row indices
    fread(ja, sizeof(int), (nrows + 1), file_bin); // column indices
    fread(a, sizeof(double), nnz, file_bin); // values

    //for(int i =0 ;i < nnz; i++) std::cout << ia[i] << ", "; std::cout << std::endl; std::cout << std::endl;
    //for(int i =0 ;i < (nrows + 1); i++) std::cout << ja[i] << ", "; std::cout << std::endl; std::cout << std::endl;

    int ncol = nrows; // Assuming square matrix
    int len = ceil(log10(0.1 + nnz + 1)) + 1;
    int nperline = std::min((int)floor(80.0/len), ncol+1);
    int ptr_len = len;
    int ptr_nperline = nperline;
    int ptrcrd = floor(ncol / nperline) + 1;
    char ptrfmt[100];
    sprintf(ptrfmt, "(%dI%d)", nperline, len);
    
    // Compute row index format
    nperline = std::min((int)floor(80.0/len), nnz);
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
        int nrow = nrows; // Assuming nrow to be nrows, adjust as necessary
        rhscrd = floor((nrhs * nrow - 1) / nperline) + 1;
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
    fprintf(fid, "%14i%14i%14i%14i\n", nrows, ncol, nnz, nrhs); // Assuming nnzeros is nnz and nrow is nrows

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
            fprintf(fid, format_string.c_str(), ja[j]);
        }
        fprintf(fid, "\n");
    }

    // Preparing the format string
    std::string formatStr = "%" + std::to_string(ind_len) + "d";

    // Loop to print all elements with ind_nperline elements per line
    for (int i = 0; i < nnz; ++i) {
        if (i > 0 && i % ind_nperline == 0) {
            fprintf(fid, "\n");
        }
        fprintf(fid, formatStr.c_str(), ia[i]);
    }

    // Add a newline if the last line isn't full
    if (nnz % nperline != 0) {
        fprintf(fid, "\n");
    }

    // Check if a newline is needed at the end
    if (nnz==ind_nperline) {
        fprintf(fid, "\n");
    }

    // Checking the initialization of the 'a' array (Add this before your code to print the values)
    //for (int i = 0; i < nnz; ++i) {
    //    std::cout << "a[" << i << "] = " << a[i] << std::endl;
    //}

    // numerical values of nonzero elements of the matrix
    if (job >= 2) {
        for (int i = 0; i < nnz; ++i) {
            fprintf(fid, c_valfmt.c_str(), a[i]);
            if ((i + 1) % nperline == 0) {
                fprintf(fid, "\n");
            }
        }

        // Add a newline if the last line isn't full
        if (nnz % nperline != 0) {
            fprintf(fid, "\n");
        }
    }
    
    fclose(fid);
    fclose(file_bin);
    free(a);
    free(ia);
    free(ja);

    return 0;
}



/*
    free(rows);
    free(columns);
    free(values);
    free(ia_csc);
    free(ja_csc);
    fclose(file);

    return 0;
}


*/


    //for(int i=0; i <=nrows; i++) std::cout << ja[i] << ", "; std::cout << std::endl; //csc formatted
    //for(int i=0; i < nnz; i++) std::cout << ia[i] << ", "; std::cout << std::endl; //row index!
 

void csc_format_to_csr_format(int nrows, int nnz, int *rows, int *columns, int *ia, int *ja){


}

void csr_format_to_csc_format(int nrows, int nnz, int *rows, int *columns, int *ia, int *ja){


}

int sort_columns(int *rows, int *columns, int n, int *sortedColumns, double *values) {


    // Step 1: Create a vector of indices [0, 1, 2, ...]
    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) {
        indices[i] = i;
    }
    
    // Step 2: Sort the indices vector based on the values in the rows array
    std::sort(indices.begin(), indices.end(), [&rows](int a, int b) {
        return rows[a] < rows[b];
    });

    // Step 3: Create new arrays to hold the sorted rows and columns
    int* sortedRows = new int[n];
    

    // Allocate a new array to hold the copy
    double* copy_values = new double[n];

    // Copy the values from the original array to the new array
    std::memcpy(copy_values, values, n * sizeof(double));

    // Step 4: Use the sorted indices to rearrange rows and columns
    for (int i = 0; i < n; ++i) {
        sortedRows[i] = rows[indices[i]];
        sortedColumns[i] = columns[indices[i]];
        values[i] = copy_values[indices[i]];
    }
    free(copy_values);

/*
     std::cout << " " << std::endl;

    // Step 5: Print the sorted arrays
    for (int i = 0; i < n; ++i) {
        std::cout << sortedRows[i] << ' ';
    }
    std::cout << '\n';

    for (int i = 0; i < n; ++i) {
        std::cout << sortedColumns[i] << ' ';
    }
    std::cout << '\n';
*/
    return 0;
}

void csc_lower_triangular_column_format(int nrows, int nnz, int *rows, int *columns, int *ia, int *ja){

    int k=1, mod = columns[0];
    ja[0] = 1;
    for(int i=0; i <nnz; i++){
        if(mod==columns[i]){
            ja[k] = ja[k] + 1;
        }else{
            mod = columns[i];
            ja[k] += ja[k-1];
            k++;
            ja[k] = ja[k] + 1;
        }
    }
    ja[nrows] = ja[nrows-1]+1;
    for(int i=0; i <nnz; i++) ia[i] = rows[i]+1;

}

void csr_lower_triangular_row_format(int nrows, int nnz, int *rows, int *columns, int *ia, int *ja, double *values) {

    int* sortedColumns = new int[nnz];
    sort_columns(rows, columns, nnz, sortedColumns, values);

    int k=1, mod = rows[0];
    ia[0] = 1;
    for(int i=0; i <nnz; i++){
        if(mod==rows[i]){
            ia[k] = ia[k] + 1;
        }else{
            mod = rows[i];
            ia[k] += ia[k-1];
            k++;
            ia[k] = ia[k] + 1;
        }

    }
    ia[nrows] = ia[nrows-1]+1;
    for(int i=0; i <nnz; i++) ja[i] = sortedColumns[i]+1;
}
