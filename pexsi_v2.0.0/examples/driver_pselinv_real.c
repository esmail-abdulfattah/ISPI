
/**
 * @file driver_pselinv_real.c
 * @brief Example for using the driver interface for parallel selected
 * inversion of a real symmetric matrix.
 *
 *
 * @date 2013-11-10 Original version.
 * @date 2014-01-26 Change the interface.
 * @date 2014-04-01 Compatible with the interface at version 0.7.0.
 * @date 2016-09-10 Compatible with the interface at version 0.10.0
 */
#include  <stdio.h>
#include  <stdlib.h>
#include <stdbool.h>

#include  "c_pexsi_interface.h"
//#include  <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



// Function to convert CSR upper triangular to full matrix
void upperToFull(int nrows, double *values, int *col_indices, int *row_ptr, 
                 double **full_values, int **full_row_ptr, int **full_col_indices) {


    //Print the first 20 elements
    printf("\nFirst 100 elements:\n");
    for (int i = 0; i < 50 && i < 100; i++) {
        printf("Value: %lf, Col index: %d, Row ptr: %d\n", values[i], col_indices[i], row_ptr[i]);
    }

    int nnz = (row_ptr[nrows]-1)*2 - nrows; // Non-zero elements in upper triangular
    printf("nnz: %d\n", nnz);

    // Allocating space for the full matrix, assuming the worst case where the matrix is dense
    *full_values = (double *)malloc(nnz * sizeof(double));
    *full_row_ptr = (int *)malloc((nrows + 1) * sizeof(int));
    *full_col_indices = (int *)malloc(nnz * sizeof(int));

    int full_nnz = 0; // Non-zero count for full matrix
    (*full_row_ptr)[0] = 1;

    for (int i = 0; i < nrows; i++) {
        for (int j = row_ptr[i]; j < row_ptr[i+1]; j++) {
            int col = col_indices[j];

            // Copying the upper triangular part
            (*full_values)[full_nnz] = values[j - row_ptr[i]];
            (*full_col_indices)[full_nnz] = col;
            full_nnz++;

            // If it's not a diagonal element, add the symmetrical lower triangular part
            if (i != col) {
                (*full_values)[full_nnz] = values[j - row_ptr[i]];
                (*full_col_indices)[full_nnz] = i;
                full_nnz++;
            }
        }
        (*full_row_ptr)[i+1] = full_nnz + 1; // Adjust the row pointers
    }

    // Printing the first few elements
    printf("\nFirst 100 elements:\n");
    for (int i = 0; i < 100 && i < full_nnz; i++) {
        printf("Value: %lf, Col index: %d, Row ptr: %d\n", (*full_values)[i], (*full_col_indices)[i], (*full_row_ptr)[i]);
    }

    // Print the last 10 elements
   // printf("\nLast 10 elements:\n");
   // for (int i = full_nnz - 10; i < full_nnz; i++) {
   //     printf("Value: %lf, Col index: %d, Row ptr: %d\n", (*full_values)[i], (*full_col_indices)[i], (*full_row_ptr)[i]);
   // }
  

}

void read_file(const char *filename, int *nrows, int *nnz, int **ia, int **ja, double **a, double **x, double **b) {

    FILE* file_bin = fopen(filename, "rb");

    if (file_bin == NULL) {
        perror("Error opening the file");
        return;
    }
 
    fread(nrows, sizeof(int), 1, file_bin);
    fread(nnz, sizeof(int), 1, file_bin);

    // Allocate memory based on the read values
    *a = (double*)malloc(*nnz * sizeof(double));
    *ia = (int*)malloc(*nnz * sizeof(int)); // row pointers
    *ja = (int*)malloc((*nrows + 1) * sizeof(int)); // column indices

    printf("%d %d\n", *nrows, *nnz);
    
    *x = (double*)malloc(*nrows * sizeof(double));
    *b = (double*)malloc(*nrows * sizeof(double));

    fread(*ia, sizeof(int), *nnz, file_bin); // column indices
    fread(*ja, sizeof(int), *nrows + 1, file_bin); // column pointers
    fread(*a, sizeof(double), *nnz, file_bin); // values
    fread(*x, sizeof(double), *nrows, file_bin); // values
    fread(*b, sizeof(double), *nrows, file_bin); // values

    fclose(file_bin);
}

void read_csr_matrix_using_bin(const char *matrix_path, int *nrows, int *nnz, int **ia, int **ja, double **a) {
    double *x, *b;  // Since x and b are not used in subsequent calls, they are local to this function.
    read_file(matrix_path, nrows, nnz, ia, ja, a, &x, &b);
    
    // Free x and b since they're not being used outside this function.
    free(x);
    free(b);
}

int default_interface(int argc, char **argv) 
{

  //upcxx-run -n 1 ./driver_pselinv_real_serverSYMPACK /home/abdulfe/donttouch/pexsi_sympack_last_test/pexsi_v2.0.0/examples/test/Q_pexsi1_num_0_dim_3779.matrix

  int mpirank, mpisize;
  int nrows = 0;     // initialize
  int nnz = 0;       // initialize
  int *colptr = NULL; // initialize 
  

  int           nnzLocal;                     
  int           numColLocal;                  
  int*          colptrLocal;                  
  int*          rowindLocal;                  
  double*       AnzvalLocal;
  double*       AinvnzvalLocal;
  int           nprow, npcol;
  int           info;
  char*         Rfile;   

  int           i, j, irow, jcol;
  int           numColLocalFirst, firstCol;

  double tsymfact1, tsymfact2, tinvert1, tinvert2;

  #ifdef WITH_SYMPACK
  symPACK_Init(&argc, &argv);
  #else
      fprintf(stderr, "Error: symPACK support not enabled.\n");
      exit(1);
  #endif

  MPI_Comm_rank( MPI_COMM_WORLD, &mpirank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpisize );

  #if defined(PROFILE) || defined(PMPI)
    TAU_PROFILE_INIT(argc, argv);
    TAU_PROFILE_SET_CONTEXT(MPI_COMM_WORLD);
  #endif



  /* Below is the data used for the toy g20 matrix */

  npcol               = 1;
  nprow               = mpisize;
  Rfile               = "lap2dr.matrix";
  //Rfile               = "laplace128full.csr";

  if (argc > 1)
      Rfile = argv[1];

  // Read the matrix 
  ReadDistSparseMatrixFormattedHeadInterface(
      Rfile,
      &nrows,
      &nnz,
      &nnzLocal,
      &numColLocal,
      MPI_COMM_WORLD );

  

  if( mpirank == 0 ){
    printf("On processor 0...\n");
    printf("nrows       = %d\n", nrows );
    printf("nnz         = %d\n", nnz );
    printf("nnzLocal    = %d\n", nnzLocal );
    printf("numColLocal = %d\n", numColLocal );


    int nrows = 0;
    int nnz = 0;
    int *ptr = NULL;
    int *ind = NULL;
    double *val = NULL;
    double *full_val = NULL;
    int *full_ptr = NULL;
    int *full_ind = NULL;

    //const char* matrix_path = "/home/abdulfe/donttouch/pexsi_sympack_last_test/pexsi_v2.0.0/examples/test/Q_pardiso_csr_num_1_dim_7479.bin";
    // read_csr_matrix_using_bin(matrix_path, &nrows, &nnz, &ptr, &ind, &val);
    
    // Print the first 20 elements
    //printf("\nFirst 20 elements:\n");
    //for (int i = 0; i < 20 && i < 100; i++) {
    //    printf("Value: %lf, Col index: %d, Row ptr: %d\n", val[i], ind[i], ptr[i]);
    //}


   // upperToFull(nrows, val, ptr, ind, &full_val, &full_ptr, &full_ind);

  }



  // Allocate memory 
  colptrLocal = (int*)malloc( (numColLocal+1) * sizeof(int) );
  rowindLocal = (int*)malloc( nnzLocal * sizeof(int) );
  AnzvalLocal = (double*)malloc( nnzLocal * sizeof(double) );
  AinvnzvalLocal = (double*)malloc( nnzLocal * sizeof(double) );
   
  // Read the matrix 
  ReadDistSparseMatrixFormattedInterface(
    Rfile,
    nrows,
    nnz,
    nnzLocal,
    numColLocal,
    colptrLocal,
    rowindLocal,
    AnzvalLocal,
    MPI_COMM_WORLD );

  if (mpirank == 0) printf("\nInitialize PEXSI\n");
  // Initialize PEXSI 

  PPEXSIOptions  options;
  INLA_Options inla_options;

  PPEXSISetDefaultOptions( &options );
  options.npSymbFact = 1;
  options.solver = 0;
  #ifdef WITH_SYMPACK
  options.solver = 1;
  #endif
  options.symmetricStorage = 1;
  options.ordering = 4;
  options.verbosity = 2;

  //options.isSymbolicFactorize = 1;
  /*
    case 0: "PARMETIS" //not working
    case 1: "METIS"
    case 2: "MMD"
    case 3: "NATURAL"
    case 4: "AMD"
    case 5: "PTSCOTCH" //not working
    case 6: "SCOTCH"
  */

  PPEXSIPlan   plan;

  plan = PPEXSIPlanInitialize( 
      MPI_COMM_WORLD, 
      nprow,
      npcol,
      mpirank, 
      &info );

  if (mpirank == 0) printf("\nPPEXSILoadRealHSMatrix\n");

  PPEXSILoadRealHSMatrix( 
      plan, 
      options,
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      AnzvalLocal,
      1,     // S is identity
      NULL,  // S is identity
      &info );

  if (mpirank == 0) tsymfact1 = MPI_Wtime(); //tsymfact1 = omp_get_wtime();
  
  if (mpirank == 0) printf("\nPEXSISymbolicFactorizeRealSymmetricMatrix\n");
      
    PPEXSISymbolicFactorizeRealSymmetricMatrix( plan,options,&info );



    if (mpirank == 0)
          printf("\n DONE \n");

    if (mpirank == 0) tsymfact2 = MPI_Wtime(); //tsymfact2 = omp_get_wtime();
      
    if (mpirank == 0) tinvert1 = MPI_Wtime(); //tinvert1 = omp_get_wtime();


    PPEXSISelInvRealSymmetricMatrix (
        plan,
        options,
        AnzvalLocal,
        AinvnzvalLocal,
        &info );
        

    if (mpirank == 0) tinvert2 = MPI_Wtime(); //tinvert2 = omp_get_wtime();
          

    if( info != 0 ){
      if( mpirank == 0 ){
        printf("PSelInv routine gives info = %d. Exit now.\n", info );
        printf("The error message is in logPEXSI* files.\n" );
      }
      MPI_Finalize();
      return info;
    }


    // The first processor output the diagonal elements in natural order
    
    if( mpirank == 0 ){
      numColLocalFirst = nrows / mpisize;
      firstCol         = mpirank * numColLocalFirst;
      for( j = 0; j < numColLocal; j++ ){
        jcol = firstCol + j + 1;
        for( i = colptrLocal[j]-1; 
            i < colptrLocal[j+1]-1; i++ ){
          irow = rowindLocal[i];
          if( irow == jcol ){
            printf("Ainv[%5d,%5d] = %15.10e\n", 
                irow, irow,
                AinvnzvalLocal[i]);
          }
        }
      } // for (j)
    }

    // Clean up 

    PPEXSIPlanFinalize( 
        plan,
        &info );
    
    if( mpirank == 0 ){ 
      printf("\nAll calculation is finished. Exit the program.\n");
    }

    if (mpirank == 0)
    {
      printf("Time needed for symb. fact. : %e \n", tsymfact2-tsymfact1);
      printf("Time needed for sel. invers.: %e \n", tinvert2-tinvert1);
    }


  // Deallocate memory 
  free( colptrLocal );
  free( rowindLocal );
  free( AnzvalLocal );
  free( AinvnzvalLocal );

  
  #ifdef WITH_SYMPACK
  symPACK_Finalize();
  #else
  MPI_Finalize();
  #endif

  return 0;
}

typedef struct {

    //
    //-------------------> MPI World
    //
    MPI_Comm comm;
    int mpirank, mpisize;
    int npcol;              
    int nprow;             
    //----------------------------------------------


    //
    //-------------------> Formats and Data:
    //
    int dim;
    
    int nnzLocal, numColLocal;  // Initialized variables
    int *colptrLocal, *rowindLocal, *numRead_vector;
    double *AnzvalLocal;

    double *AinvnzvalLocal;
    //----------------------------------------------

    //
    //-------------------> Phases for repetitions
    //
    bool phase_0, phase_1, phase_2, phase_3, phase_4, phase_5, phase_6, phase_7, phase_8;
    //----------------------------------------------


    //
    //-------------------> PEXSI Objects
    //
    int nnz_pexsi, *pexsi_csc_p_j, *pexsi_csc_i, *pexsi_x_order, num_calls_pexsi;
    double *pexsi_x;
    PPEXSIOptions  pexsi_options;
    PPEXSIPlan   plan;
    int info;
    //----------------------------------------------

    //
    //-------------------> PEXSI Objects
    //
    int nnz_sympack, *sympack_csc_p_j, *sympack_csc_i, *sympack_x_order, num_calls_sympack;
    double *sympack_x;
    bool use_sympack2D;
    //----------------------------------------------

    bool opt_symbolicfactorization, opt_factorization, opt_logdeterminant, opt_solve_Ax_b, opt_solve_Ly_z, opt_solve_LTy_b, opt_inverse;
    double* rhs_b;
    double* sol_x;
    int num_rhs;

} INLA_OBJECT;

void default_inla_object(INLA_OBJECT *obj, MPI_Comm comm, int nnz, int dim);
int get_formats(int argc, char **argv);
void inla_phase1(INLA_OBJECT *inla_options);
void inla_phase2(INLA_OBJECT *inla_options);
void inla_phase5(INLA_OBJECT *inla_options);
void inla_phase6(INLA_OBJECT *inla_options);
char* replace_substring(const char *str, const char *old, const char *newStr);

int inla_sympack_pexsi_interface(INLA_OBJECT *inla_options, MPI_Comm comm){

  //if comm is changed in inla_options, then call this phase 1 again.
  inla_phase1(inla_options);
  if(inla_options->mpirank==0) printf("\nPhase 1 is completed\n\n");

    if( inla_options->mpirank == 0 ){
    printf("On processor 0...\n");
    printf("nrows       = %d\n", inla_options->dim );
    printf("nnz         = %d\n", inla_options->nnz_pexsi );
    printf("nnzLocal    = %d\n", inla_options->nnzLocal );
    printf("numColLocal = %d\n", inla_options->numColLocal );
  }

  //memory allocation for the values.
  inla_phase2(inla_options); 
  if(inla_options->mpirank==0) printf("\nPhase 2 is completed\n");

  //goal is to call this once
  //INLA_ReadDistSparseMatrixFormattedInterface_phase_1(
  inla_phase3(
      inla_options->dim,
      inla_options->nnz_pexsi, 
      inla_options->nnzLocal,
      inla_options->numColLocal, 
      inla_options->colptrLocal, 
      inla_options->rowindLocal, 
      inla_options->comm ,
      inla_options->pexsi_csc_p_j,
      inla_options->pexsi_csc_i,
      inla_options->numRead_vector,
      &inla_options->phase_3);
  if(inla_options->mpirank==0) printf("\nPhase 3 is completed\n");

  //goal is to call this everytime you change the values.
  //INLA_ReadDistSparseMatrixFormattedInterface_phase_2(
  inla_phase4(
     inla_options->comm, 
     inla_options->pexsi_x,
     inla_options->AnzvalLocal,
     inla_options->nnzLocal,
     inla_options->numRead_vector,
     &inla_options->phase_4);
  if(inla_options->mpirank==0) printf("\nPhase 4 is completed\n");

  //call this if mpi world changes.
  inla_phase5(inla_options);
  if(inla_options->mpirank==0) printf("\nPhase 5 is completed\n");

  //goal is to call this everytime you change the values.
  inla_phase6(inla_options); //PPEXSILoadRealHSMatrix
  if(inla_options->mpirank==0) printf("\nPhase 6 is completed\n\n");



/*
  //this should be called once.
  INLA_PPEXSISymbolicFactorizeRealSymmetricMatrix( inla_options->plan, inla_options->pexsi_options,&inla_options->info);

  //this should be called many times because values in the matrix changes.
  INLA_PPEXSIFactorizeRealSymmetricMatrix(inla_options->plan, inla_options->pexsi_options, inla_options->AnzvalLocal, 
                                          inla_options->AinvnzvalLocal, &inla_options->info);

  //Here are the solve steps and logdeterminant computation.
  INLA_PPEXSI_Solve_LogDeterminant_RealSymmetricMatrix(inla_options->plan, inla_options->pexsi_options, inla_options->AnzvalLocal, 
                                          inla_options->AinvnzvalLocal, &inla_options->info,
                                          inla_options->sol_x, inla_options->rhs_b, inla_options->num_rhs,
                                          inla_options->opt_logdeterminant, inla_options->opt_solve_Ax_b, 
                                          inla_options->opt_solve_Ly_z, inla_options->opt_solve_LTy_b, inla_options->dim);

  //Only called when the partial inversion is needed.
  INLA_PPEXSISelInvRealSymmetricMatrix(inla_options->plan, inla_options->pexsi_options, inla_options->AnzvalLocal, inla_options->AinvnzvalLocal, &inla_options->info );
*/                


    PPEXSISymbolicFactorizeRealSymmetricMatrix( inla_options->plan, inla_options->pexsi_options, &inla_options->info);


    PPEXSISelInvRealSymmetricMatrix (inla_options->plan, inla_options->pexsi_options, inla_options->AnzvalLocal, inla_options->AinvnzvalLocal,  &inla_options->info);

  if( inla_options->info != 0 ){
    if(inla_options->mpirank == 0 ){
      printf("PSelInv routine gives info = %d. Exit now.\n", inla_options->info );
      printf("The error message is in logPEXSI* files.\n" );
    }
    MPI_Finalize();
    return inla_options->info;
  }

    // The first processor output the diagonal elements in natural order
    
    int           i, j, irow, jcol;
    int           numColLocalFirst, firstCol;

    if(inla_options->mpirank == 0 ){
      numColLocalFirst = inla_options->dim / inla_options->mpisize;
      firstCol         =inla_options->mpirank * numColLocalFirst;
      for( j = 0; j < inla_options->numColLocal; j++ ){
        jcol = firstCol + j + 1;
        for( i = inla_options->colptrLocal[j]-1; 
            i < inla_options->colptrLocal[j+1]-1; i++ ){
          irow = inla_options->rowindLocal[i];
          if( irow == jcol ){
            if(irow==1 || irow==inla_options->dim) printf("Ainv[%5d,%5d] = %15.10e\n", irow, irow, inla_options->AinvnzvalLocal[i]);
          }
        }
      } // for (j)
    }

    // Clean up 
    PPEXSIPlanFinalize( inla_options->plan, &inla_options->info );

  return 1;

}

int new_interface(int argc, char **argv);

int main(int argc, char **argv) {

  bool onlysympack = false;

  if(!onlysympack){

    if (argc < 2) {
            fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
            return 1;
        }

        const char *filename = argv[1];
        const char *dot_bin = ".bin";
        size_t len_filename = strlen(filename);
        size_t len_bin = strlen(dot_bin);

        if ((len_filename > len_bin && strcmp(filename + len_filename - len_bin, dot_bin) == 0)){

            //upcxx-run -n 1 ./driver_pselinv_real_serverSYMPACK /home/abdulfe/donttouch/pexsi_sympack_last_test/pexsi_v2.0.0/examples/matrices/group1/Q_pardiso_csr_num_101_dim_8763.bin
            new_interface(argc, argv);     // Use upper CSR if the file ends with .bin
        
        } else {

            //upcxx-run -n 1 ./driver_pselinv_real_serverSYMPACK /home/abdulfe/donttouch/pexsi_sympack_last_test/pexsi_v2.0.0/examples/matrices/group1/Q_pexsi_num_1_dim_8763.matrix
            default_interface(argc, argv); // Use full CSC formats otherwise: the default pexsi interface
        
        }

  }else{
    
    //the default sympack interface: this we add just for comparison.
    //use for example: upcxx-run -n 1 ./driver_pselinv_real_serverSYMPACK /home/abdulfe/donttouch/pexsi_sympack_last_test/pexsi_v2.0.0/examples/matrices/group1/Q_sympack_csc_num_1_dim_8763.rb
    sympack_interface(argc, argv); 
  
  }
   

    return 0;
}



void default_inla_object(INLA_OBJECT *obj, MPI_Comm comm, int nnz, int dim) {
    
    //-------------------------------------------
    //Format Phases:
    //-------------------------------------------

    obj->dim = dim;
    obj->nnz_pexsi = 2*nnz - dim;
    obj->nnz_sympack = nnz;

    obj->pexsi_csc_p_j = NULL;
    obj->pexsi_csc_i = NULL;
    obj->pexsi_x = (double*)malloc((2*nnz-dim) * sizeof(double));
    obj->pexsi_x_order = NULL;
    obj->num_calls_pexsi = 0;

    obj->sympack_csc_p_j = NULL;
    obj->sympack_csc_i = NULL;
    obj->sympack_x = (double*)malloc(nnz * sizeof(double));
    obj->sympack_x_order = NULL;
    obj->use_sympack2D = false; //TODO?! :( PEXSI doesn't support it)

    obj->rhs_b = NULL;
    obj->sol_x = NULL;

    obj->num_calls_sympack = 0;

    obj->phase_0 = false; // or true, depending on your desired default
    obj->phase_1 = false;
    obj->phase_2 = false;
    obj->phase_3 = false;
    obj->phase_4 = false;
    obj->phase_5 = false;
    obj->phase_6 = false;

    obj->comm = comm;
    MPI_Comm_rank( obj->comm, &obj->mpirank );
    MPI_Comm_size( obj->comm, &obj->mpisize );
    obj->nnzLocal = 0;
    obj->numColLocal = 0;
    obj->numRead_vector = (int*)malloc(obj->mpisize * sizeof(int));
    obj->use_sympack2D = false;

    //-------------------------------------------
    //PEXSI OPTIONS:
    //-------------------------------------------
    PPEXSISetDefaultOptions(&obj->pexsi_options);
    obj->pexsi_options.npSymbFact = 1;
    obj->pexsi_options.solver = 0;

    #ifdef WITH_SYMPACK
    obj->pexsi_options.solver = 1;
    #endif

    obj->pexsi_options.symmetricStorage = 1;
    obj->pexsi_options.ordering = 4; //{'0':"PARMETIS", '1':"METIS", '2':"MMD", '3':"NATURAL", '4':"AMD", '5':"PTSCOTCH", '6':"SCOTCH"}
    obj->pexsi_options.verbosity = 1;
    obj->pexsi_options.isSymbolicFactorize = 1;

    obj->npcol = 1;   //npcol*nprow = mpisize, factorization depends only on npcol?
    obj->nprow = obj->mpisize;             

    obj->opt_symbolicfactorization = true;
    obj->opt_logdeterminant = true;
    obj->opt_solve_Ax_b = true;
    obj->opt_solve_Ly_z = true; //L: lower traingle we get from cholesky
    obj->opt_solve_LTy_b = true; //LT: upper traingle we get from cholesky
    obj->opt_inverse = true;

    obj->num_rhs = 1; //number of right hand sides.

}

int new_interface(int argc, char **argv) {

  #ifdef WITH_SYMPACK
  symPACK_Init(&argc, &argv);
  #else
      fprintf(stderr, "Error: symPACK support not enabled.\n");
      exit(1);
  #endif

  MPI_Comm comm = MPI_COMM_WORLD; // Corrected this line

  double* a = NULL;
  int* ia = NULL;
  int* p_ja = NULL;
  int nnz;
  int dim;

  double* rhs_b = NULL;
  double* sol_x = NULL;
  int num_rhs;

  INLA_OBJECT inla_options;
  inla_options.phase_0 = false;

  //Replace "101" with "202" in argv[1]
  char *filename1 = strdup(argv[1]);
  char *filename2 = strdup(replace_substring(argv[1], "101", "202"));


  get_upper_csr_bin_c_filename(argv[1], &a, &ia, &p_ja, &nnz, &dim, true, true, &rhs_b, &sol_x, &num_rhs);
 
      //
      //----------------> Create inla objects with some default values
      //
      default_inla_object(&inla_options, comm, nnz, dim);
      //here you can create function to change default parameters. 

      //
      //----------------> Get sympack format from the upper csr format - should be called once - sympack_x is not changed.
      //
      inla_format_interface_sympack(&a, &ia, &p_ja, &nnz, &dim, 
                        &inla_options.sympack_csc_i, &inla_options.sympack_csc_p_j, &inla_options.sympack_x_order, &inla_options.sympack_x,
                        &inla_options.num_calls_sympack);

      //
      //----------------> Get pexsi format from the upper csr format - should be called once - pexsi_x is not changed.
      //
      inla_format_interface_pexsi(&a, &ia, &p_ja, &nnz, &dim, 
                          &inla_options.pexsi_csc_i, &inla_options.pexsi_csc_p_j, &inla_options.pexsi_x_order, &inla_options.pexsi_x,
                          &inla_options.num_calls_pexsi);

      
      printf("Phase 0 is completed");

      //
      //----------------> update x for sympack format in case we are using sympack: this step is to be called everytime the values x change. Since the pattern is the same - 
      //
      inla_format_interface_sympack(&a, &ia, &p_ja, &nnz, &dim, 
                          &inla_options.sympack_csc_i, &inla_options.sympack_csc_p_j, &inla_options.sympack_x_order, &inla_options.sympack_x,
                          &inla_options.num_calls_sympack);
                          
      //
      //----------------> update x for pexsi format in case we are using pexsi: this step is to be called everytime the values x change. Since the pattern is the same - 
      //
      inla_format_interface_pexsi(&a, &ia, &p_ja, &nnz, &dim, 
                        &inla_options.pexsi_csc_i, &inla_options.pexsi_csc_p_j, &inla_options.pexsi_x_order, &inla_options.pexsi_x,
                        &inla_options.num_calls_pexsi, &inla_options.nnz_pexsi, &inla_options.dim);
  

      //
      //----------------> get right hand side and "solution x" for comparison.
      //
      if(inla_options.opt_solve_Ax_b || inla_options.opt_solve_LTy_b || inla_options.opt_solve_Ly_z){

        inla_options.num_rhs = num_rhs;
        inla_options.rhs_b = (double*)malloc(inla_options.dim * inla_options.num_rhs* sizeof(double));
        inla_options.sol_x = (double*)malloc(inla_options.dim * inla_options.num_rhs* sizeof(double));
        inla_options.rhs_b = rhs_b;
        inla_options.sol_x = sol_x;
      }

      inla_sympack_pexsi_interface(&inla_options, comm);
  

  
  // Deallocate memory 
  free( inla_options.colptrLocal );
  free( inla_options.rowindLocal );
  free( inla_options.AnzvalLocal );
  free( inla_options.AinvnzvalLocal );

  #ifdef WITH_SYMPACK
  symPACK_Finalize();
  #else
  MPI_Finalize();
  #endif

  return 1;
}

void inla_phase1(INLA_OBJECT *inla_options){

  //----------------------------------------------------------
  //phase 1: get nnzLocal and numColLocal (no values for now)
  //----------------------------------------------------------
  
  inla_phase1_call(&inla_options->dim, &inla_options->nnz_pexsi, inla_options->pexsi_csc_p_j, &inla_options->nnzLocal, &inla_options->numColLocal, inla_options->comm);
  
  // Allocate memory 
  inla_options->colptrLocal = (int*)malloc((inla_options->numColLocal + 1) * sizeof(int));
  inla_options->rowindLocal = (int*)malloc(inla_options->nnzLocal * sizeof(int));
  inla_options->AnzvalLocal = (double*)malloc(inla_options->nnzLocal * sizeof(double));
  
  inla_options->phase_1 = true;
}

void inla_phase2(INLA_OBJECT *inla_options) {

    if(inla_options->AinvnzvalLocal == NULL){
      inla_options->AinvnzvalLocal = (double*)malloc(inla_options->nnzLocal * sizeof(double));
    }

    inla_options->phase_2 = true;

}

void inla_phase5(INLA_OBJECT *inla_options) {

  inla_options->plan = PPEXSIPlanInitialize( inla_options->comm, inla_options->nprow,inla_options->npcol, inla_options->mpirank, &inla_options->info );
  inla_options->phase_5 = true;
}

void inla_phase6(INLA_OBJECT *inla_options) {

  PPEXSILoadRealHSMatrix( 
      inla_options->plan, 
      inla_options->pexsi_options,
      inla_options->dim,
      inla_options->nnz_pexsi,
      inla_options->nnzLocal,
      inla_options->numColLocal,
      inla_options->colptrLocal,
      inla_options->rowindLocal,
      inla_options->AnzvalLocal,
      1,     // S is identity
      NULL,  // S is identity
      &inla_options->info );

  inla_options->phase_6 = true;

}

char* replace_substring(const char *str, const char *old, const char *newStr) {
    
    char *result;
    int i, cnt = 0;
    int newlen = strlen(newStr);
    int oldlen = strlen(old);

    // Counting the number of times old word occurs in the string
    for (i = 0; str[i] != '\0'; i++) {
        if (strstr(&str[i], old) == &str[i]) {
            cnt++;
            i += oldlen - 1;
        }
    }

    // Making new string of enough length
    result = (char *)malloc(i + cnt * (newlen - oldlen) + 1);

    i = 0;
    while (*str) {
        // Compare the substring with the result
        if (strstr(str, old) == str) {
            strcpy(&result[i], newStr);
            i += newlen;
            str += oldlen;
        }
        else
            result[i++] = *str++;
    }

    result[i] = '\0';
    return result;
}







/*




int inla_functions() 
{

  int nnzLocal;
  int numColLocal;                  


          
  int*          colptrLocal;                  
  int*          rowindLocal;                  
  double*       AnzvalLocal;
  double*       AinvnzvalLocal;
  int           nprow, npcol;
  int           info;
  char*         Rfile;   

  int           i, j, irow, jcol;
  int           numColLocalFirst, firstCol;

  double tsymfact1, tsymfact2, tinvert1, tinvert2;







  // Below is the data used for the toy g20 matrix 

  npcol               = 1;
  nprow               = mpisize;
  Rfile               = "lap2dr.matrix";
  //Rfile               = "laplace128full.csr";

  if (argc > 1)
      Rfile = argv[1];

  // Read the matrix 
  if(true){

      ReadDistSparseMatrixFormattedHeadInterface(
          Rfile,
          &nrows,
          &nnz,
          &nnzLocal,
          &numColLocal,
          MPI_COMM_WORLD );

  }else{
      
      //TODO: the input should be only csr format for inla
      inla_phase1_call(&nrows, &nnz, colptr, &nnzLocal, &numColLocal, MPI_COMM_WORLD);
      
  }

  if( mpirank == 0 ){
    printf("On processor 0...\n");
    printf("nrows       = %d\n", nrows );
    printf("nnz         = %d\n", nnz );
    printf("nnzLocal    = %d\n", nnzLocal );
    printf("numColLocal = %d\n", numColLocal );
  }

  // Allocate memory 
  colptrLocal = (int*)malloc( (numColLocal+1) * sizeof(int) );
  rowindLocal = (int*)malloc( nnzLocal * sizeof(int) );
  AnzvalLocal = (double*)malloc( nnzLocal * sizeof(double) );
  AinvnzvalLocal = (double*)malloc( nnzLocal * sizeof(double) );

  if(true){
    
    // Read the matrix 
    ReadDistSparseMatrixFormattedInterface(
      Rfile,
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      AnzvalLocal,
      MPI_COMM_WORLD );

  }else{
    
    //TODO: the input should be only csr format for inla
    // Read the matrix 
    INLA_ReadDistSparseMatrixFormattedInterface(
        Rfile,
        nrows,
        nnz,
        nnzLocal,
        numColLocal,
        colptrLocal,
        rowindLocal,
        AnzvalLocal,
        MPI_COMM_WORLD );
  }

  if (mpirank == 0) printf("\nInitialize PEXSI\n");
  // Initialize PEXSI 

  PPEXSIOptions  options;
  INLA_Options inla_options;

  PPEXSISetDefaultOptions( &options );
  options.npSymbFact = 1;
  options.solver = 0;
  #ifdef WITH_SYMPACK
  options.solver = 1;
  #endif
  options.symmetricStorage = 1;
  options.ordering = 1;
  options.verbosity = 2;

    //options.isSymbolicFactorize = 1;
    //
    //  case 0: "PARMETIS" //not working
    //  case 1: "METIS"
    //  case 2: "MMD"
    //  case 3: "NATURAL"
    //  case 4: "AMD"
    //  case 5: "PTSCOTCH" //not working
    //  case 6: "SCOTCH"
    //

  PPEXSIPlan   plan;

  plan = PPEXSIPlanInitialize( 
      MPI_COMM_WORLD, 
      nprow,
      npcol,
      mpirank, 
      &info );

  if (mpirank == 0) printf("\nPPEXSILoadRealHSMatrix\n");

  PPEXSILoadRealHSMatrix( 
      plan, 
      options,
      nrows,
      nnz,
      nnzLocal,
      numColLocal,
      colptrLocal,
      rowindLocal,
      AnzvalLocal,
      1,     // S is identity
      NULL,  // S is identity
      &info );

  if (mpirank == 0) tsymfact1 = MPI_Wtime(); //tsymfact1 = omp_get_wtime();
  
  if (mpirank == 0) printf("\nPEXSISymbolicFactorizeRealSymmetricMatrix\n");

  if(false){
      
    PPEXSISymbolicFactorizeRealSymmetricMatrix( plan,options,&info );



    if (mpirank == 0)
          printf("\n DONE \n");

    if (mpirank == 0) tsymfact2 = MPI_Wtime(); //tsymfact2 = omp_get_wtime();
      
    if (mpirank == 0) tinvert1 = MPI_Wtime(); //tinvert1 = omp_get_wtime();


    PPEXSISelInvRealSymmetricMatrix (
        plan,
        options,
        AnzvalLocal,
        AinvnzvalLocal,
        &info );
        

    if (mpirank == 0) tinvert2 = MPI_Wtime(); //tinvert2 = omp_get_wtime();
          

    if( info != 0 ){
      if( mpirank == 0 ){
        printf("PSelInv routine gives info = %d. Exit now.\n", info );
        printf("The error message is in logPEXSI* files.\n" );
      }
      MPI_Finalize();
      return info;
    }


    // The first processor output the diagonal elements in natural order
    
    if( mpirank == 0 ){
      numColLocalFirst = nrows / mpisize;
      firstCol         = mpirank * numColLocalFirst;
      for( j = 0; j < numColLocal; j++ ){
        jcol = firstCol + j + 1;
        for( i = colptrLocal[j]-1; 
            i < colptrLocal[j+1]-1; i++ ){
          irow = rowindLocal[i];
          if( irow == jcol ){
            printf("Ainv[%5d,%5d] = %15.10e\n", 
                irow, irow,
                AinvnzvalLocal[i]);
          }
        }
      } // for (j)
    }

    // Clean up 

    PPEXSIPlanFinalize( 
        plan,
        &info );
    
    if( mpirank == 0 ){ 
      printf("\nAll calculation is finished. Exit the program.\n");
    }

    if (mpirank == 0)
    {
      printf("Time needed for symb. fact. : %e \n", tsymfact2-tsymfact1);
      printf("Time needed for sel. invers.: %e \n", tinvert2-tinvert1);
    }

  }else{

      INLA_PHASES_INTERFACE( 
          plan,
          options,
          &info,
          AnzvalLocal,
          AinvnzvalLocal);

  }

  // Deallocate memory 
  free( colptrLocal );
  free( rowindLocal );
  free( AnzvalLocal );
  free( AinvnzvalLocal );

  
  #ifdef WITH_SYMPACK
  symPACK_Finalize();
  #else
  MPI_Finalize();
  #endif

  return 0;
}
*/


void read_some_matrix(int *dim, int *nnz, int **colptr, int argc, char **argv){

    char* Rfile = argv[1];

    FILE* file = fopen(Rfile, "r");
    if (!file) {
        perror("Error opening file");
        return;
    }

    // Read matrix dimensions and nnz
    int zero_check;
    fscanf(file, "%d %d %d %d", dim, dim, nnz, &zero_check);

    // Allocate memory for column pointers and row indices
    *colptr = (int*)malloc((*dim + 1) * sizeof(int));
    int *rowind = (int*)malloc(*nnz * sizeof(int));
    if (!(*colptr) || !rowind) {
        perror("Error allocating memory");
        fclose(file);
        return;
    }

    // Read column pointers
    for (int i = 0; i <= *dim; i++) {
        fscanf(file, "%d", &((*colptr)[i]));
    }

    // Read row indices
    for (int i = 0; i < *nnz; i++) {
        fscanf(file, "%d", &rowind[i]);
    }

    fclose(file);

  /*
      void inla_phaseold0(INLA_OBJECT *inla_options, int argc, char **argv) {
      
      read_some_matrix(&inla_options->dim, &inla_options->nnz_pexsi, &inla_options->colptr, argc, argv);
      inla_options->phase_0 = true;
  }*/

}
