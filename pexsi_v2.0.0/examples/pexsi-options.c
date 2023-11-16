    double        spin;   //NumSpin, default is 2.0
    double        temperature;  //Temperature, in the same unit as H 
    double        gap; //Spectral gap. **Note** This can be set to be 0 in most cases.
    double        deltaE; //@brief  An upper bound for the spectral radius of \f$S^{-1} H\f$.
    double        muMin0; //@brief  Initial guess of lower bound for mu.
    double        muMax0; //@brief  Initial guess of upper bound for mu.
    double        mu0; //@brief  Initial guess for mu (for the solver) (AG)
    double        muInertiaTolerance; //@brief  Stopping criterion in terms of the chemical potential for the inertia counting procedure.
    double        muInertiaExpansion; //@brief  If the chemical potential is not in the initial interval, the interval is expanded by muInertiaExpansion.
    double        muPEXSISafeGuard; //@brief  Safe guard criterion in terms of the chemical potential to reinvoke the inertia counting procedure.
    double        numElectronPEXSITolerance; //@brief  Stopping criterion of the %PEXSI iteration in terms of the number of electrons compared to numElectronExact.

    int           numPole; //@brief  Number of terms in the pole expansion.
    int           isInertiaCount; //@brief  Whether inertia counting is used at the very beginning.
    int           maxPEXSIIter; //Maximum number of %PEXSI iterations after each inertia counting procedure.
    int           matrixType; //@brief  matrixType (global) Type of input H and S matrices. {'0': Real symmetric (default), '1': General complex matrices (not implemented yet)}
    int           isSymbolicFactorize; //@brief  Whether to perform symbolic factorization.
    int           isConstructCommPattern; // @brief  Whether to construct PSelInv communication pattern.
    int           solver; //@brief  Solver used to do the factorization prior to the selected inversion. {'0': SuperLU_DIST, '1': symPACK (For symmetric matrices only)}
    int           symmetricStorage; //@brief  Storage space used by the Selected Inversion algorithm for symmetric matrices. {'0'   : Non symmetric storage, '1'   : Symmetric storage (lower memory usage)}


int           solver = 1; 
int           symmetricStorage = 1;
int           ordering;
int           rowOrdering;
int           npSymbFact = 1;
int           symmetric = 0;
int           verbosity;

    /** 
     * @brief  Ordering strategy for factorization and selected
     * inversion. When SuperLU is used:  
     * - = 0   : Parallel ordering using ParMETIS/PT-SCOTCH (PARMETIS
     *   option in SuperLU_DIST).
     * - = 1   : Sequential ordering using METIS (METIS_AT_PLUS_A
     *   option in SuperLU_DIST).
     * - = 2   : Multiple minimum degree ordering (MMD_AT_PLUS_A
     *   option in SuperLU_DIST).
     * When symPACK is used:
     * - = 0   : Parallel ordering using PT-SCOTCH.
     * - = 1   : Sequential ordering using SCOTCH.
     * - = 2   : Multiple minimum degree ordering.
     * - = 3   : Approximate minimum degree ordering.
     * - = 4   : Parallel ordering using PARMETIS.
     * - = 5   : Sequential ordering using METIS.
     */ 
    int           ordering;
    /** 
     * @brief  row permutation strategy for factorization and selected
     * inversion.  
     * - = 0   : No row permutation (NOROWPERM
     *   option in SuperLU_DIST).
     * - = 1   : Make diagonal entry larger than off diagonal ( LargeDiag
     *   option in SuperLU_DIST).
     */ 
    int           rowOrdering;
    /** 
     * @brief  Number of processors for PARMETIS/PT-SCOTCH.  Only used
     * if the ordering == 0.
     */ 
    int           npSymbFact;
    /** 
     * @brief  Matrix structure.
     * - = 0   : Unsymmetric matrix
     * - = 1   : Symmetric matrix (default).
     */ 
    int           symmetric;
    /** 
     * @brief  Transpose.
     * - = 0   : Factor non transposed matrix (default).
     * - = 1   : Factor transposed matrix.
     */ 
    int           transpose;
    /** 
     * @brief  The pole expansion method to be used.
     * - = 1   : Cauchy Contour Integral method used.
     * - = 2   : Moussa optimized method.
     */ 
    int           method;
    /** 
     * @brief  The point parallelizaion of PEXSI.
     * - = 2  : Recommend two points parallelization
     */ 
    int           nPoints;
    /** 
     * @brief  The driver version of the PEXSI.
     * - = 1.0.0 : Latest version is 1.0.0 
     */ 
    //char         driverVersion[10] ;//= "1.0.0";
    /** 
     * @brief  The level of output information.
     * - = 0   : No output.
     * - = 1   : Basic output (default)
     * - = 2   : Detailed output.
     */ 
    int           verbosity;

    /** 
     * @brief  The error message of the PEXSI.
     * - = 0   : No error.
     * - = 1   : Smaller mu has bigger number of electrons
     */ 
     int         iFLAG;
