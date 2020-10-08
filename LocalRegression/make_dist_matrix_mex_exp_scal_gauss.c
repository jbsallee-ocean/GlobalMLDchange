
#include "mex.h"
#include <math.h>
#include <omp.h>



void mexFunction(int nlhs, mxArray *plhs[], 
    int nrhs, const mxArray *prhs[])
{
 int jm, jmk, j, k, m, mkj;
 double *data0, *data1, *data2;
 double *scal;
 double tmp;

   
   
    /* Find the dimensions of the data */
    m = mxGetM(prhs[0]);
 

    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix(m, m, mxREAL);


    /* Retrieve the input data */
    data0 = mxGetPr(prhs[0]);
    data1 = mxGetPr(prhs[1]);
    scal = mxGetPr(prhs[2]);
    

    /* Create a pointer to the output data */
    data2 = mxGetPr(plhs[0]);


     /* Put data in the output array */
    for (j = 0; j < m; j++)
    {
        jm = j*m;
    #pragma omp for        
    for (k = 0; k < j; k++) 
    {
    jmk = jm + k;
    mkj = m*k+j;
    data2[jmk] = data1[j] - data1[k];
    data2[jmk] = data2[jmk] * data2[jmk];
    data2[mkj] = data0[j] - data0[k];
    data2[mkj] = data2[mkj] * data2[mkj]; 
    data2[jmk] = data2[jmk] + data2[mkj];
    data2[jmk] = data2[jmk] / scal[0];
    data2[mkj] = data2[jmk];
    }
    }
   
   
}

