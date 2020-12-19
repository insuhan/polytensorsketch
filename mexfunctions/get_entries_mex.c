#include "mex.h"

double get_entries(double *x, int i, int j, mwSize nrows, mwSize ncols)
{
    mwSize k;
    double out;
    
    out = 0.0;
    for (k = 0; k < ncols; k++){
        out += x[i + k*nrows] * x[j + k*nrows];        
    }
    return out;    
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *x;
    double *idx1;
    double *idx2;

    mwSize ncols;
    mwSize nrows;
    mwSize nidx;
    mwSize r1,r2,c1,c2;
    mwSize i,j;
    
    double *y;
    
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:get_entries:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:get_entries:nlhs","One output required.");
    }
    x = mxGetPr(prhs[0]);
    idx1 = mxGetPr(prhs[1]);
    idx2 = mxGetPr(prhs[2]);
    
    nrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    
    r1 = mxGetM(prhs[1]);
    c1 = mxGetN(prhs[1]);
    
    r2 = mxGetM(prhs[2]);
    c2 = mxGetN(prhs[2]);
    
    if(c1!=c2) {
        mexErrMsgIdAndTxt("MyToolbox:get_entries:nlhs","Index size should be same.");
    }
    if(r1!=1 || r2!=1) {
        mexErrMsgIdAndTxt("MyToolbox:get_entries:nlhs","Index size should be same.");
    }
    
//     printf("nrows = %d, ncols = %d\n", nrows, ncols);
//     printf("rc1 = %d, %d, | rc2 = %d, %d\n", r1,c1,r2,c2);
    
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)c1,mxREAL);
    y = mxGetPr(plhs[0]);
    
    for (i=0; i<c1; i++)
    {
        y[i] = get_entries(x, (int)(idx1[i]) - 1, (int)(idx2[i]) - 1, nrows, ncols);
    }    
}


