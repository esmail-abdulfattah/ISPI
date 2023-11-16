
get_full <- function(sparse_Q){
  
  # Extract the upper triangular part of the sparse matrix
  upperTri <- triu(sparse_Q)
  
  # Extract the lower triangular part of the sparse matrix
  lowerTri <- tril(sparse_Q)
  
  # Combine the upper and lower triangular parts to form a complete sparse matrix
  fullTriSparse <- upperTri + t(upperTri) + lowerTri + t(lowerTri)
  diag(fullTriSparse) <- diag(sparse_Q)
  
  # Print the non-zero elements of the complete sparse matrix
  return(fullTriSparse)
}
write_lower_mat_csc <- function(Q, filename, num, export_b) {
  
  Q <- t(Q)
  
  dim_Q = dim(Q)[1]
  b = numeric(dim_Q)
  x = round(rnorm(dim_Q, 3, sd = 2),3)
  b = get_full(Q)%*%x
  
  ia <- Q@i + 1
  ja <- Q@p + 1
  values <- Q@x
  
  # Print the first 10 elements of ia array
  cat("ia[1:10]:\n")
  cat(ia[1:10], sep = ", ", fill = TRUE)
  cat("\n\n")
  
  # Print the first 10 elements of ja array
  cat("ja[1:10]:\n")
  cat(ja[1:10], sep = ", ", fill = TRUE)
  cat("\n\n")
  
  # Print the first 10 elements of values array
  cat("values[1:10]:\n")
  cat(values[1:10], sep = ", ", fill = TRUE)
  cat("\n")
  
  
  nnz <- length(values)
  dim <- dim(Q)[1]
  
  # Open the binary file for writing
  filename <- paste0(filename, "_num_", as.integer(num), "_dim_", dim, ".bin")
  print(filename)
  fwrite <- file(filename, "wb")
  
  # Write header information
  writeBin(object = as.integer(dim), con = fwrite, size = 4)
  # Print the function name
  cat("write_lower_mat_csc - R - Function\n")
  
  # Print the nnz value
  cat("nnz: ", as.integer(nnz), "\n")
  
  # Print the dim value
  cat("dim: ", as.integer(dim), "\n")
  
  
  writeBin(object = as.integer(nnz), con = fwrite, size = 4)
  
  # Write matrix data in CSC format
  writeBin(object = as.integer(ia), con = fwrite, size = 4)
  writeBin(object = as.integer(ja), con = fwrite, size = 4)
  writeBin(object = as.double(values), con = fwrite, size = 8)
  writeBin(object = as.double(x), con = fwrite, size = 8)
  writeBin(object = as.double(b), con = fwrite, size = 8)
  
  #print(as.double(x))
  #print(as.double(b))
  print(paste("first element", Q[1,1], " - last element", Q[dim,dim]))
  # Close the file
  close(con = fwrite)
}
export_to_pexsi2 <- function(matrix, filename, num) {
  
  # Extract the indices and values
  rows <- matrix@i
  cols <- matrix@p
  values <- matrix@x
  
  # Extract the dimensions
  dims <- matrix@Dim
  
  file_path <- paste0(filename, "_num_", as.integer(num), "_dim_", dims[1], ".matrix")
  
  # Open a connection to write to the file
  con <- file(file_path, "wt")
  
  # Write the dimensions, number of non-zero elements, and 0
  num_non_zero <- length(values)
  cat(paste(dims[1], dims[2], num_non_zero, 0), "\n", file = con)
  
  # Write the 'p' vector (cols)
  cat(paste0(cols + 1, collapse = " "), "\n", file = con)
  
  # Write the 'i' vector (rows)
  cat(paste0(rows + 1, collapse = " "), "\n", file = con)
  
  # Write the values, each indented by 8 spaces
  sapply(values, function(val) {
    cat(sprintf("        %f", val), "\n", file = con)
  })
  
  # Close the connection
  close(con)
}
# Function to export dense matrix to mat.matrix format
export_to_pexsi1 <- function(M, filename, num) {
  # Initialize variables
  rowind <- c()
  colind <- c()
  values <- c()
  rows <- nrow(M)
  cols <- ncol(M)

  filename <- paste0(filename, "_num_", as.integer(num), "_dim_", rows, ".matrix")
  
  # Loop through the dense matrix to get the row indices, column indices and non-zero values
  for (j in 1:cols) {
    for (i in 1:rows) {
      if (M[i, j] != 0) {
        rowind <- c(rowind, i)
        colind <- c(colind, j)
        values <- c(values, M[i, j])
      }
    }
  }
  
  # Number of non-zero values
  nnz <- length(values)
  
  # Create column pointers
  colptr <- c(1, cumsum(table(colind)) + 1)
  
  # Open a file connection
  con <- file(filename, "w")
  
  # Write the first line (rows, columns, nnz, 0)
  write(paste(rows, cols, nnz, 0), file = con)
  
  # Write the second line (colptr)
  write(paste(colptr, collapse = " "), file = con)
  
  # Write the third line (rowind)
  write(paste(rowind, collapse = " "), file = con)
  
  # Convert the values to double precision format
  double_values <- sprintf("%.16e", values)
  
  # Write the remaining lines (values)
  write(paste(double_values, collapse = "\n"), file = con)
  
  # Close the file connection
  close(con)
}

# export_to_pexsi2 <- function(M, filename) {
#   # Initialize variables
#   rowind <- c()
#   colind <- c()
#   values <- c()
#   rows <- nrow(M)
#   cols <- ncol(M)
#   
#   # Loop through the dense matrix to get the row indices, column indices, and non-zero values
#   # Only consider upper diagonal elements including the diagonal itself
#   for (j in 1:cols) {
#     for (i in 1:rows) {
#       if (i <= j && M[i, j] != 0) {
#         rowind <- c(rowind, i)
#         colind <- c(colind, j)
#         values <- c(values, M[i, j])
#       }
#     }
#   }
#   
#   # Number of non-zero values
#   nnz <- length(values)
#   
#   # Create column pointers
#   colptr <- c(1, cumsum(table(colind)) + 1)
#   
#   # Open a file connection
#   con <- file(filename, "w")
#   
#   # Write the first line (rows, columns, nnz, 0)
#   write(paste(rows, cols, nnz, 0), file = con)
#   
#   # Write the second line (colptr)
#   write(paste(colptr, collapse = " "), file = con)
#   
#   # Write the third line (rowind)
#   write(paste(rowind, collapse = " "), file = con)
#   
#   # Convert the values to double precision format
#   double_values <- sprintf("%.16e", values)
#   
#   # Write the remaining lines (values)
#   write(paste(double_values, collapse = "\n"), file = con)
#   
#   # Close the file connection
#   close(con)
# }
