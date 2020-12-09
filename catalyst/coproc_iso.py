
#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt=0
forceOutputAtFirstCall=False

# Global screenshot output options
imageFileNamePadding=0
rescale_lookuptable=False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays=False

# a root directory under which all Catalyst output goes
rootDirectory=''

# makes a cinema D index table
make_cinema_table=False

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.8.1
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.8.1

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.8.1
      #
      # To ensure correct image size when batch processing, please search 
      # for and uncomment the line `# renderView*.ViewSize = [*,*]`

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # get the material library
      materialLibrary1 = GetMaterialLibrary()

      LoadPalette(paletteName='GradientBackground')

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1556, 1044]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [3.0, 1.5, 0.5]
      renderView1.StereoType = 'Crystal Eyes'
      renderView1.CameraPosition = [-4.233567488134634, -7.966379736173562, 6.40855107227724]
      renderView1.CameraFocalPoint = [3.0000000000000013, 1.5000000000000002, 0.4999999999999999]
      renderView1.CameraViewUp = [0.24374792047590432, 0.37261896006415096, 0.8954005036096322]
      renderView1.CameraFocalDisk = 1.0
      renderView1.CameraParallelScale = 3.4418828593064434
      renderView1.BackEnd = 'OSPRay raycaster'
      renderView1.OSPRayMaterialLibrary = materialLibrary1

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='isosurfaces/isosurface_%t.png', freq=1, fittoscreen=0, magnification=1, width=1000, height=670, cinema={}, compression=5)
      renderView1.ViewTime = datadescription.GetTime()

      SetActiveView(None)

      # ----------------------------------------------------------------
      # setup view layouts
      # ----------------------------------------------------------------

      # create new layout object 'Layout #1'
      layout1 = CreateLayout(name='Layout #1')
      layout1.AssignView(0, renderView1)

      # ----------------------------------------------------------------
      # restore active view
      SetActiveView(renderView1)
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a producer from a simulation input
      input = coprocessor.CreateProducer(datadescription, 'input')

      ## create a new 'Contour'
      #contour1 = Contour(Input=input)
      #contour1.ContourBy = ['POINTS', 'Q']
      #contour1.Isosurfaces = [1.0]
      #contour1.PointMergeMethod = 'Uniform Binning'
      # create a new 'VTKm Contour'
      contour1 = VTKmContour(Input=input)
      contour1.ContourBy = ['POINTS', 'Q']
      contour1.Isosurfaces = [1.0]
      contour1.ComputeScalars = 1

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from input
      inputDisplay = Show(input, renderView1, 'GeometryRepresentation')

      # trace defaults for the display properties.
      inputDisplay.Representation = 'Outline'
      inputDisplay.ColorArrayName = [None, '']
      inputDisplay.OSPRayScaleArray = 'Pressure'
      inputDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      inputDisplay.SelectOrientationVectors = 'None'
      inputDisplay.ScaleFactor = 0.609375
      inputDisplay.SelectScaleArray = 'None'
      inputDisplay.GlyphType = 'Arrow'
      inputDisplay.GlyphTableIndexArray = 'None'
      inputDisplay.GaussianRadius = 0.03046875
      inputDisplay.SetScaleArray = ['POINTS', 'Pressure']
      inputDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      inputDisplay.OpacityArray = ['POINTS', 'Pressure']
      inputDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      inputDisplay.DataAxesGrid = 'GridAxesRepresentation'
      inputDisplay.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
      inputDisplay.OSPRayScaleFunction.Points = [-0.014005885448460252, 0.0, 0.5, 0.0, 1.1572820062885005, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      inputDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      inputDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.1757813367477812e-38, 1.0, 0.5, 0.0]

      # show data from contour1
      contour1Display = Show(contour1, renderView1, 'GeometryRepresentation')

      # get color transfer function/color map for 'U'
      uLUT = GetColorTransferFunction('U')
      uLUT.RGBPoints = [0.21536746402820187, 0.0, 0.0, 0.0, 0.7290007786732375, 0.901960784314, 0.0, 0.0, 1.2426340933182733, 0.901960784314, 0.901960784314, 0.0, 1.499450750640791, 1.0, 1.0, 1.0]
      uLUT.ColorSpace = 'RGB'
      uLUT.NanColor = [0.0, 0.498039215686, 1.0]
      uLUT.ScalarRangeInitialized = 1.0

      # trace defaults for the display properties.
      contour1Display.Representation = 'Surface'
      contour1Display.ColorArrayName = ['POINTS', 'U']
      contour1Display.LookupTable = uLUT
      contour1Display.OSPRayScaleArray = 'Q'
      contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      contour1Display.SelectOrientationVectors = 'None'
      contour1Display.ScaleFactor = 0.12835979461669922
      contour1Display.SelectScaleArray = 'Q'
      contour1Display.GlyphType = 'Arrow'
      contour1Display.GlyphTableIndexArray = 'Q'
      contour1Display.GaussianRadius = 0.006417989730834961
      contour1Display.SetScaleArray = ['POINTS', 'Q']
      contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
      contour1Display.OpacityArray = ['POINTS', 'Q']
      contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
      contour1Display.DataAxesGrid = 'GridAxesRepresentation'
      contour1Display.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
      contour1Display.OSPRayScaleFunction.Points = [-0.014005885448460252, 0.0, 0.5, 0.0, 1.1572820062885005, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      contour1Display.ScaleTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      contour1Display.OpacityTransferFunction.Points = [1.0, 0.0, 0.5, 0.0, 1.000244140625, 1.0, 0.5, 0.0]

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for uLUT in view renderView1
      uLUTColorBar = GetScalarBar(uLUT, renderView1)
      uLUTColorBar.Title = 'U'
      uLUTColorBar.ComponentTitle = ''

      # set color bar visibility
      uLUTColorBar.Visibility = 1

      # show color legend
      contour1Display.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get opacity transfer function/opacity map for 'U'
      uPWF = GetOpacityTransferFunction('U')
      uPWF.Points = [0.21536746402820187, 0.0, 0.5, 0.0, 1.499450750640791, 1.0, 0.5, 0.0]
      uPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(contour1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  if requestSpecificArrays:
    arrays = [['Pressure', 0], ['Q', 0], ['U', 0], ['V', 0], ['W', 0]]
    coprocessor.SetRequestedArrays('input', arrays)
  coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,forceOutputAtFirstCall)

  if rootDirectory:
      coprocessor.SetRootDirectory(rootDirectory)

  if make_cinema_table:
      coprocessor.EnableCinemaDTable()

  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
