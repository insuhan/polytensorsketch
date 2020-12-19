# Polynomial Tensor Sketch
MATLAB implementation for [Polynomial Tensor Sketch for Element-wise Function of Low-Rank Matrix](https://arxiv.org/abs/1905.11616) (ICML 2020)

## Installation
Run ```install.m``` to compile mex files in ```./mexfunctions``` and add paths:
```matlab
>> install
```

## Usage
Run ```run_kernel_approx.m``` with a proper input for kernel approximation tasks:
```matlab
>> run_kernel_approx('synthetic')
segment dataset is loaded, n=2310, d=19
PTS (coreset) error: 0.001349
RFF           error: 0.022135
>> run_kernel_approx('segment') 
segment dataset is loaded, n=2310, d=19
PTS (coreset) error: 0.000176
RFF           error: 0.003045
>> run_kernel_approx('satimage') 
satimage dataset is loaded, n=4435, d=36
PTS (coreset) error: 0.000097
RFF           error: 0.008614
```