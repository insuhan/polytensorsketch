clc
mex -largeArrayDims -lmwblas -DUSE_BLAS CountSketchMex.c
mex -largeArrayDims get_entries_mex.c