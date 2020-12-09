#include <iostream>
#include <algorithm>

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkMPI.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkRectilinearGrid.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>

#include <cuda_runtime.h>


extern "C"
{
namespace
{
struct GlobalVars {
  vtkCPProcessor* Processor = NULL;
  vtkMPICommunicatorOpaqueComm* comm = NULL;

#if 0
  // Note: RectilinearGrid is better for this case but vtkm contour crashes with it...
  vtkRectilinearGrid* FlowGrid;
#else
  vtkStructuredGrid* FlowGrid;
#endif
  float *xc, *yc, *zc, *xyzc;
  int nx, ny, nz;
  int procDims[2];
  int procIdx[2];
};

GlobalVars globals;

// Routine to build mesh in VTK
void BuildFlowGrid(unsigned int nx, double lx, unsigned int ny, double ly, unsigned int nz, double *zc)
{

  double dx = lx / nx;
  double dy = ly / ny;

#if 0
  // Using managed memory here to avoid some memory movement in Paraview/Catalyst. Not strictly required.
  cudaMallocManaged(&globals.xc, (nx + 2) * sizeof(float));
  cudaMallocManaged(&globals.yc, (ny + 2) * sizeof(float));
  cudaMallocManaged(&globals.zc, nz*sizeof(float));

  // Set point coordinates (offset by -1 in i,j for ghosts)
  for (int i = 0; i < nx + 2; ++i) globals.xc[i] = globals.procIdx[0] * lx + (i-1) * dx + dx / 2;
  for (int j = 0; j < ny + 2; ++j) globals.yc[j] = globals.procIdx[1] * ly + (j-1) * dy + dy / 2;
  for (int k = 0; k < nz; ++k) globals.zc[k] = zc[k];
#else

  // Using managed memory here to avoid some memory movement in Paraview/Catalyst. Not strictly required.
  cudaMallocManaged(&globals.xyzc, (nx + 2) * (ny + 2) * nz * 3 * sizeof(float));
  float* point = globals.xyzc;
  for (int k = 0; k < nz; ++k) {
    for (int j = 0; j < ny + 2; ++j) {
      for (int i = 0; i < nx + 2; ++i) {
        point[0] = globals.procIdx[0] * lx + (i-1) * dx + dx / 2;
        point[1] = globals.procIdx[1] * ly + (j-1) * dy + dy / 2;
        point[2] = zc[k];
	point += 3;
      }
    }
  }
  vtkNew<vtkFloatArray> pointArray;
  pointArray->SetNumberOfComponents(3);
  pointArray->SetArray(globals.xyzc, static_cast<vtkIdType>((nx + 2) * (ny + 2) * nz * 3), 1);
  vtkNew<vtkPoints> points;
  points->SetData(pointArray.GetPointer());
#endif

#if 0
  vtkNew<vtkFloatArray> xc, yc, zc;
  xc->SetNumberOfComponents(1);
  yc->SetNumberOfComponents(1);
  zc->SetNumberOfComponents(1);
  xc->SetArray(globals.xc, static_cast<vtkIdType>(nx + 2), 1);
  yc->SetArray(globals.yc, static_cast<vtkIdType>(ny + 2), 1);
  zc->SetArray(globals.zc, static_cast<vtkIdType>(nz), 1);

  globals.FlowGrid->SetXCoordinates(xc.GetPointer());
  globals.FlowGrid->SetYCoordinates(yc.GetPointer());
  globals.FlowGrid->SetZCoordinates(zc.GetPointer());
#else
  globals.FlowGrid->SetPoints(points);
#endif

  globals.FlowGrid->SetExtent(globals.procIdx[0] * nx, (globals.procIdx[0] + 1) * nx + 1,
		              globals.procIdx[1] * ny, (globals.procIdx[1] + 1) * ny + 1,
			      0, nz - 1);

}

// Routine to associate data to field grid
void InitializeFlowGridAttributes(unsigned int numberOfPoints,
                                  double* uData, double* vData, double *wData, double* pData, double* qcritData)
{
  // Note: Comment out any of these fields to remove from output data files
  vtkNew<vtkDoubleArray> u;
  u->SetName("U");
  u->SetNumberOfComponents(1);
  u->SetArray(uData, static_cast<vtkIdType>(numberOfPoints), 1);
  globals.FlowGrid->GetPointData()->AddArray(u.GetPointer());

  vtkNew<vtkDoubleArray> v;
  v->SetName("V");
  v->SetNumberOfComponents(1);
  v->SetArray(vData, static_cast<vtkIdType>(numberOfPoints), 1);
  globals.FlowGrid->GetPointData()->AddArray(v.GetPointer());

  vtkNew<vtkDoubleArray> w;
  w->SetName("W");
  w->SetNumberOfComponents(1);
  w->SetArray(wData, static_cast<vtkIdType>(numberOfPoints), 1);
  globals.FlowGrid->GetPointData()->AddArray(w.GetPointer());

  vtkNew<vtkDoubleArray> p;
  p->SetName("Pressure");
  p->SetNumberOfComponents(1);
  p->SetArray(pData, static_cast<vtkIdType>(numberOfPoints), 1);
  globals.FlowGrid->GetPointData()->AddArray(p.GetPointer());

  vtkNew<vtkDoubleArray> qcrit;
  qcrit->SetName("Q");
  qcrit->SetNumberOfComponents(1);
  qcrit->SetArray(qcritData, static_cast<vtkIdType>(numberOfPoints), 1);
  globals.FlowGrid->GetPointData()->AddArray(qcrit.GetPointer());

}


// Routine to build flow grid and associate velocity data arrays
void InitializeFlowGrid(unsigned int nx, double lx, unsigned int ny, double ly, unsigned int nz, double* zc,
                        double* uData, double* vData, double *wData, double* pData, double* qcritData,
                        int* procDims, int* procIdx)
{
  unsigned int numberOfPoints = (nx + 2) * (ny + 2) * nz; //includes ghost points in x, y

  for (int i = 0; i < 2; ++i) {
    globals.procDims[i] = procDims[i];
    globals.procIdx[i] = procIdx[i];
  }

  globals.nx = nx;
  globals.ny = ny;
  globals.nz = nz;

  std::cout << "Calling InitializeFlowGrid:" << std::endl;
  std::cout << "\tnumber of points (with ghost cells): " << numberOfPoints << std::endl;

#if 0
  globals.FlowGrid = vtkRectilinearGrid::New();
#else
  globals.FlowGrid = vtkStructuredGrid::New();
#endif
  BuildFlowGrid(nx, lx, ny, ly, nz, zc);
  InitializeFlowGridAttributes(numberOfPoints, uData, vData, wData, pData, qcritData);
}

void CatalystInitialize(bool active)
{
  int color = MPI_UNDEFINED;
  if (active) color = 1;

  MPI_Comm handle;
  MPI_Comm_split(MPI_COMM_WORLD, color, 0, &handle);

  if (globals.Processor == NULL and active)
  {
    globals.Processor = vtkCPProcessor::New();
    globals.comm = new vtkMPICommunicatorOpaqueComm(&handle);
    globals.Processor->Initialize(*globals.comm);

    vtkNew<vtkCPPythonScriptPipeline> pipeline;
    pipeline->Initialize("coproc.py");
    globals.Processor->AddPipeline(pipeline.GetPointer());
  }
}

void CatalystFinalize()
{
  if (globals.Processor)
  {
    globals.Processor->Delete();
    globals.Processor = NULL;
  }
  if (globals.FlowGrid)
  {
    globals.FlowGrid->Delete();
    globals.FlowGrid = NULL;
  }
}
}

// Routine to update Catalyst coprocesser pipeline.
void CatalystCoProcess(double time, unsigned int timeStep)
{

  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput("input");
  dataDescription->SetTimeData(time, timeStep);

  if (globals.Processor->RequestDataDescription(dataDescription.GetPointer()) != 0)
  {
    dataDescription->GetInputDescriptionByName("input")->SetGrid(globals.FlowGrid);
    dataDescription->GetInputDescriptionByName("input")->SetWholeExtent(0, globals.procDims[0] * globals.nx + 1,
		                                                        0, globals.procDims[1] * globals.ny + 1,
			                                                0, globals.nz - 1);
    globals.Processor->CoProcess(dataDescription.GetPointer());
  }
}
}
