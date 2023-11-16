
setwd("~/server/donttouch/pexsi_sympack_last_test/symPACK/test")

source("export_matrix.R")

library(Matrix)
set.seed(123)  # Set seed for reproducibility

n = 100
# Generate a random 9x9 matrix
mat <- matrix(rpois(n * n, lambda = 10), n, n)

# Create a positive definite matrix
positive_definite_mat <- round(t(mat) %*% mat,4)
C <- chol(positive_definite_mat)

eigen(positive_definite_mat)$values

num = dim(sparse_matrix)[1]

positive_definite_mat <- Matrix(positive_definite_mat, sparse = TRUE)
positive_definite_mat <- as(positive_definite_mat, "dgCMatrix")

#export_to_pexsi1(positive_definite_mat, "size_100a_compare1.matrix") #595.24304049573913744
export_to_pexsi1(positive_definite_mat, "Q_pexsi1", num)
export_to_pexsi2(positive_definite_mat, "Q_pexsi2", num)

#image(Qxy)
sparse_matrix <- triu(positive_definite_mat)
sparse_Q <- as(sparse_matrix, "sparseMatrix")

# Call the function with the correct name
write_lower_mat_csc(sparse_Q, "Q_sympack", num, TRUE)
system("./convert_to_rb.sh")


2*sum(log(diag(C)))
