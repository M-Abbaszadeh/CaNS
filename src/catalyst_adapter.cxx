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
  unsigned int ng[3];
};

GlobalVars globals;

// Routine to build mesh in VTK
void BuildFlowGrid(unsigned int* n, double* xc, double* yc, double* zc)
{
#if 0
  // Using managed memory here to avoid some memory movement in Paraview/Catalyst. Not strictly required.
  cudaMallocManaged(&globals.xc, n[0] * sizeof(float));
  cudaMallocManaged(&globals.yc, n[1] * sizeof(float));
  cudaMallocManaged(&globals.zc, n[2] * sizeof(float));

  for (int i = 0; i < n[0]; ++i) globals.xc[i] = xc[k];
  for (int j = 0; j < n[1]; ++j) globals.yc[j] = yc[k];
  for (int k = 0; k < n[2]; ++k) globals.zc[k] = zc[k];
#else

  // Using managed memory here to avoid some memory movement in Paraview/Catalyst. Not strictly required.
  cudaMallocManaged(&globals.xyzc, n[0] * n[1] * n[2] * 3 * sizeof(float));
  float* point = globals.xyzc;
  for (int k = 0; k < n[2]; ++k) {
    for (int j = 0; j < n[1]; ++j) {
      for (int i = 0; i < n[0]; ++i) {
        point[0] = xc[i];
        point[1] = yc[j];
        point[2] = zc[k];
	point += 3;
      }
    }
  }
  vtkNew<vtkFloatArray> pointArray;
  pointArray->SetNumberOfComponents(3);
  pointArray->SetArray(globals.xyzc, static_cast<vtkIdType>(n[0] * n[1] * n[2] * 3), 1);
  vtkNew<vtkPoints> points;
  points->SetData(pointArray.GetPointer());
#endif

#if 0
  vtkNew<vtkFloatArray> xc, yc, zc;
  xc->SetNumberOfComponents(1);
  yc->SetNumberOfComponents(1);
  zc->SetNumberOfComponents(1);
  xc->SetArray(globals.xc, static_cast<vtkIdType>(n[0]), 1);
  yc->SetArray(globals.yc, static_cast<vtkIdType>(n[1]), 1);
  zc->SetArray(globals.zc, static_cast<vtkIdType>(n[2]), 1);

  globals.FlowGrid->SetXCoordinates(xc.GetPointer());
  globals.FlowGrid->SetYCoordinates(yc.GetPointer());
  globals.FlowGrid->SetZCoordinates(zc.GetPointer());
#else
  globals.FlowGrid->SetPoints(points);
#endif

  globals.FlowGrid->SetExtent(0,n[0]-1,0,n[1]-1,0,n[2]-1);

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
void InitializeFlowGrid(unsigned int* n, unsigned int* ng, double* xc, double* yc, double* zc,
                        double* uData, double* vData, double *wData, double* pData, double* qcritData)
{
  unsigned int numberOfPoints = n[0] * n[1] * n[2]; //includes ghost points in x, y

  for (int i = 0; i < 3; ++i) {
    globals.ng[i] = ng[i];
  }

  std::cout << "Calling InitializeFlowGrid:" << std::endl;
  std::cout << "\tnumber of points (with ghost cells): " << numberOfPoints << std::endl;

#if 0
  globals.FlowGrid = vtkRectilinearGrid::New();
#else
  globals.FlowGrid = vtkStructuredGrid::New();
#endif
  BuildFlowGrid(n, xc, yc, zc);
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
    dataDescription->GetInputDescriptionByName("input")->SetWholeExtent(0, globals.ng[0] - 1,
		                                                        0, globals.ng[1] - 1,
			                                                0, globals.ng[2] - 1);
    globals.Processor->CoProcess(dataDescription.GetPointer());
  }
}
}
