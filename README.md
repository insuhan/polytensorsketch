# Polynomial Tensor Sketch
MATLAB implementation for [Polynomial Tensor Sketch for Element-wise Function of Low-Rank Matrix](https://arxiv.org/abs/1905.11616) (ICML 2020)

## Installation
Run ```install.m``` to compile mex files in ```./mexfunctions``` and add paths:
```matlab
>> install
```

## Usage
Run ```run_kernel_approx.m``` for testing the kernel approximation task:
```matlab
>> run_kernel_approx
segment dataset is loaded, n=2310, d=19
PTS (coreset) error: 0.001349
RFF           error: 0.022135
```