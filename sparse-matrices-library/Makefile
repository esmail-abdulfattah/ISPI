# Define the compiler
CXX=g++

# Compiler flags
CXXFLAGS=-std=c++11 -Wall -Wextra -fPIC

# Source files
SOURCES=sparse_matrix_formats.cpp main.cpp
# Header files (if needed for dependencies)
HEADERS=sparse_matrix_formats.hpp

# Output library name
LIBNAME=libsparseformats.so

# Shared library flags
LDFLAGS=-shared

# 'all' target
all: $(LIBNAME)

$(LIBNAME): $(SOURCES) $(HEADERS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SOURCES)

# 'clean' target
clean:
	rm -f $(LIBNAME)

# Phony targets
.PHONY: all clean

