#include <mex.h>

#include <vtkPolyData.h>
#include <vtkDoubleArray.h>

#include <MeshAnalyser.h>


void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[])
{

    double *v;
    double *f;
    int nV, nF, i, j;

    if(nrhs<2)
        mexErrMsgTxt("2 input arguments required");
        
    bool normalize = false;
    
    nV = mxGetM(prhs[0]);
    nF = mxGetM(prhs[1]);
    v = mxGetPr(prhs[0]);
    f = (double*)mxGetPr(prhs[1]);

    vtkPolyData *surf = vtkPolyData::New();
    vtkPoints *sverts = vtkPoints::New();
    vtkCellArray *sfaces = vtkCellArray::New();

    // Load the point, cell, and data attributes.
    //indexes in VTK are 0-based, in Matlab 1-based
    const int offset = -1;
    for(i=0; i<nV; i++) sverts->InsertPoint(i,v[i],v[nV+i],v[2*nV+i]);
    for(i=0; i<nF; i++)
    {
	    sfaces->InsertNextCell(3);
	    for(j=0;j<3;j++)
        {
		    int idx = (unsigned int)f[j*nF+i]+offset;
		    if(idx<0 || idx>=nV)
		        mexErrMsgTxt("wrong index");
		    sfaces->InsertCellPoint(vtkIdType(idx));
	    }
    }

    // assign the pieces to the vtkPolyData.
    surf->SetPoints(sverts);
    surf->SetPolys(sfaces);
    sverts->Delete();
    sfaces->Delete();

    MeshAnalyser* depthComputer = new MeshAnalyser(surf);
    depthComputer->ComputeTravelDepthFromClosed(normalize);
    vtkDoubleArray* depthMap = depthComputer->GetTravelDepth();

    if (nlhs > 1)
        mexErrMsgTxt("Too many output arguments");
    
    int nComp = depthMap->GetNumberOfComponents();
    int nRows = depthMap->GetNumberOfTuples();
    mxArray* m = mxCreateDoubleMatrix(nRows, 1, mxComplexity::mxREAL);
    double* pm = mxGetDoubles(m);
    
    
    for (i=0; i<nRows; i++)
    {
        double currData = 0;
        depthMap->GetTuple(i, &currData);
        pm[i] = currData;
    }


    plhs[0] = m;
}
