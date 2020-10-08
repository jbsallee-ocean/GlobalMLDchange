#include "mex.h"
#include <math.h>
#include <omp.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int jm, jmk, j, k, m, n;
    double *data0, *data1, *data2;
    double *scal, *signoise;
    
    /* Find the dimensions of the data */
    m = mxGetM(prhs[0]);
    
    
    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix(m, m, mxREAL);
    
    
    /* Retrieve the input data */
    data0 = mxGetPr(prhs[0]);
    data1 = mxGetPr(prhs[1]);
    scal  = mxGetPr(prhs[2]);
    signoise = mxGetPr(prhs[3]);
    
    /* Create a pointer to the output data */
    data2 = mxGetPr(plhs[0]);
    
    /* Compute upper triangle and copy identical values into lower triangle of matrix */
    
    for (j = 0; j < m; j++)
    {
        jm = j*m;
        #pragma omp for
                for (k = 0; k < j; k++)
                {
            jmk = jm + k;
            data2[jmk] =  data0[jmk];
	    data2[jmk] =  exp (-sqrt(data2[jmk]));
	    data2[jmk] =  data2[jmk] / sqrt( scal[j] * scal[k] );
            data2[m*k+j] = data2[jmk];
                }
    }
    
    #pragma omp for
            for (j = 0; j < m; j++)
            {
        data2[j*(m+1)] = signoise[j];
            }
    
}
