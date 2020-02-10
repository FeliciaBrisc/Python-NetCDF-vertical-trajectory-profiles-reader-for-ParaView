""" 

If you didn't already set the Python/System path for the netcdf4-python package
then set it up here 

import sys  
sys.path.append('/...') 

"""


import netCDF4
import vtk  
from netCDF4 import Dataset
import numpy as np


# The paraview.util.vtkAlgorithm module provides VTKPythonAlgorithmBase, the base class
# for all python-based vtkAlgorithm subclasses in VTK and decorators used to
# 'register' the algorithm with ParaView along with information about the GUI.
from paraview.util.vtkAlgorithm import *


"""
Code by Felicia Brisc (CEN University of Hamburg), distributed under a BSD 3-Clause License

A Python filter for ParaView (www.paraview.org). Reads vertical trajectory files in NetCDF format and the corresponding time evolution

Version 1.0 

The reader requires the external module netcdf4-python

The examples made available by Kitware at the link below have provided the starting point for this reader
https://gitlab.kitware.com/paraview/paraview/blob/master/Examples/Plugins/PythonAlgorithm/PythonAlgorithmExamples.py

"""


def createModifiedCallback(anobject):
    import weakref
    weakref_obj = weakref.ref(anobject)
    anobject = None
    def _markmodified(*args, **kwars):
        o = weakref_obj()
        if o is not None:
            o.Modified()
    return _markmodified



# Decorators used to add the reader 
# @smproxy.source(name="PythonNetCDFProfilesReader", label="Python-based NetCDF Profiles Reader")
# @smhint.xml("""<ReaderFactory extensions="nc" file_description="Numpy NetCDF Profiles files" />""")
@smproxy.reader(name="PythonNetCDFProfilesReader", label="Python-based NetCDF Profiles Reader",
                extensions="nc", file_description="NetCDF Profiles files")



class PythonNetCDFProfilesReader(VTKPythonAlgorithmBase):
    """A reader that reads a NetCDF vertical profile file. The NetCDF file needs to  a "time" dimension, so 
    the data is treated as a temporal dataset"""
    
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkPolyData')
        self.filename = None
        self.nc_dataset = None
        
        self.newPoints = vtk.vtkPoints()
        
        #this will be the polyline output, drawn on screen
        self.pdo = vtk.vtkPolyData()
        self.polyIndex = 0
        
              
        #the mandatory dimensions - will be put separately from the other (generic) variables  
        self.lat = None
        self.lon = None
        self.height = None
        self.time = None
        
        #self.nc_variables = [] hmmm
        
        self.nc_point_variables = []
        self.nc_cell_variables  = []
        
        self.nc_dimensions = []
        self.timesteps = None #that previous stuff 
  	
        self.timeStepsCount = 0
        #will be possibly used in the future to figure out wether the timesteps are played forward or backward
        self.timeCompare = 0 
        self.heightCount = 0
        self.numPoints = 0

        from vtkmodules.vtkCommonCore import vtkDataArraySelection
        self.arrayselection = vtkDataArraySelection()
        self.arrayselection.AddObserver("ModifiedEvent", createModifiedCallback(self))
 
 
        
        
    def build_nc_array(self, nc_var):
        vtk_array = None
   
        if nc_var.dtype.type is np.float64:
            vtk_array = vtk.vtkDoubleArray()
        if nc_var.dtype.type is np.float32:
            vtk_array = vtk.vtkFloatArray()
        if nc_var.dtype.type is np.int8 or nc_var.dtype.type is np.int16:
            vtk_array = vtk.vtkShortArray()             
        if nc_var.dtype.type is np.int32 or nc_var.dtype.type is np.int64:
            vtk_array = vtk.vtkIntArray()
        
        #extract actual data from the data set         
        nc_var_data = nc_var[:]
        
        vtk_array.SetName(nc_var.name)
        vtk_array.SetNumberOfComponents(1)
        vtk_array.SetNumberOfTuples(nc_var_data.size)
        vtk_array.SetArray(nc_var_data.data, nc_var_data.size, True)
        vtk_array.array = nc_var_data.data
      
        return vtk_array
  
        

    def get_nc_data(self, requested_time=None):
        
        if self.nc_dataset is not None:
            #if requested_time is not None:                
                #return self.nc_dataset[self.nc_dataset["time"]==requested_time]             
            
            
            return self.nc_dataset
  
        if self.filename is None:
            # Note, exceptions are totally fine!
            raise RuntimeError("No filename specified")

        #read in the data
        self.nc_dataset = Dataset(self.filename, 'r+')
      
        #find and load the mandatory dimensions (time, lon, lat, height) and the generic variables of the dataset
        nc_info = [var for var in self.nc_dataset.variables.keys()]
        nc_info = np.append(nc_info, [dim for dim in self.nc_dataset.dimensions.keys()])
      
        for i in nc_info: 
            nc_attribs = self.nc_dataset.variables[i].ncattrs()  
          
            standard_match = ['standard_name']
            matching = [s for s in nc_attribs if any(xs in s for xs in standard_match)]
                       
            if len(matching)>0 and matching[0] == 'standard_name':            
                if self.nc_dataset.variables[i].getncattr(matching[0]) == 'time':
                    self.time = self.nc_dataset.variables[i][:]
                elif self.nc_dataset.variables[i].getncattr(matching[0]) == 'latitude':
                    self.lat = self.nc_dataset.variables[i][:]
                elif self.nc_dataset.variables[i].getncattr(matching[0]) == 'longitude':
                    self.lon = self.nc_dataset.variables[i][:]
                elif self.nc_dataset.variables[i].getncattr(matching[0]) == 'height':
                    self.height = self.nc_dataset.variables[i][:]
                else: 
                    #check if it's point data of cell data (i.e. point data will depend on time and height, cell data only on time)
                    #point data
                    if self.nc_dataset.variables[i].ndim == 2: 
                        self.nc_point_variables.append(self.build_nc_array(self.nc_dataset.variables[i]))
                    
                    #cell data    
                    else:  
                        self.nc_cell_variables.append(self.build_nc_array(self.nc_dataset.variables[i]))
            
            #and if there is any variable that doesn't have the standard_name specified
            else:
                #check if it's point data of cell data (i.e. point data will depend on time and height, cell data only on time)
                #point data
                if self.nc_dataset.variables[i].ndim == 2:
                    self.nc_point_variables.append(self.build_nc_array(self.nc_dataset.variables[i]))
                #cell data
                else:  
                    self.nc_cell_variables.append(self.build_nc_array(self.nc_dataset.variables[i]))  

        
        # check if we have all four mandatory dimensions
        assert self.time is not None, "The time dimension is missing! The data will not be displayed correctly."
        assert self.lat is not None, "The latitude dimension is missing! The data will not be displayed correctly."
        assert self.lon is not None, "The longitude dimension is missing! The data will not be displayed correctly."
        assert self.height is not None, "The height dimension is missing! The data will not be displayed correctly."
 
        self.timeStepsCount = len(self.time.data)
        self.heightCount = len(self.height.data)
        self.numPoints = self.timeStepsCount * self.heightCount
        
        print("timeSteps: " , self.timeStepsCount)
        print("heightCount:" , self.heightCount)
        print("numPoints: ", self.numPoints) 
        

        self.timesteps = None
 
        if len(self.time) > 0:
            self.timesteps = self.time
            
        for i in self.nc_dataset.variables.keys(): 
            
        	self.arrayselection.AddArray(i)	

        return self.get_nc_data(requested_time)

    def get_timesteps(self):
        self.get_nc_data()
        return self.timesteps.tolist() if self.timesteps is not None else None

    def get_update_time(self, outInfo):
        executive = self.GetExecutive()
        timesteps = self.get_timesteps()
        if timesteps is None or len(timesteps) == 0:
            return None
        elif outInfo.Has(executive.UPDATE_TIME_STEP()) and len(timesteps) > 0:
            utime = outInfo.Get(executive.UPDATE_TIME_STEP())
            dtime = timesteps[0]
            for atime in timesteps:
                if atime > utime:
                    return dtime
                else:
                    dtime = atime
            return dtime
        else:
            assert(len(timesteps) > 0)
            return timesteps[0]

    def _get_array_selection(self):
        return self.arrayselection
       
        
        
    def drawPolyLine(self, currentAnimationTime):    

        aPolyLine = vtk.vtkPolyLine()
        aPolyLine.GetPointIds().SetNumberOfIds(self.heightCount)
        
        for j in range(0, self.heightCount):
            x = self.lon[currentAnimationTime]
            y = self.lat[currentAnimationTime]
            z = self.height[j]
            self.newPoints.InsertPoint(self.polyIndex+j, x,y,z) 
            aPolyLine.GetPointIds().SetId(j, self.polyIndex+j) 
            
        self.polyIndex += self.heightCount
             
        #update the output polyline
        self.pdo.SetPoints(self.newPoints)        
        self.pdo.InsertNextCell(aPolyLine.GetCellType(), aPolyLine.GetPointIds())        
        self.timeCompare = currentAnimationTime


        
    #GUI
    @smproperty.stringvector(name="FileName")
    @smdomain.filelist()
    @smhint.filechooser(extensions="nc", file_description="NetCDF trajectory profiles files")
    def SetFileName(self, name):
        """Specify filename for the file to read."""
        if self.filename != name:
            self.filename = name
            self.nc_dataset = None
            self.timesteps = None
            self.Modified()

    @smproperty.doublevector(name="TimestepValues", information_only="1", si_class="vtkSITimeStepsProperty")
    def GetTimestepValues(self):
        return self.get_timesteps()

    # The VTK array selection API allows users to choose which arrays to
    # load. smproperty.dataarrayselection() exposes this in ParaView
    # This method *must* return a `vtkDataArraySelection` instance
    @smproperty.dataarrayselection(name="Arrays")
    def GetDataArraySelection(self):
        return self._get_array_selection()



    #RequestInformation is called for the initial loading of the data
    def RequestInformation(self, request, inInfoVec, outInfoVec):

        executive = self.GetExecutive()
        outInfo = outInfoVec.GetInformationObject(0)
        outInfo.Remove(executive.TIME_STEPS())
        outInfo.Remove(executive.TIME_RANGE())
        
        #get data information and do initial data loading
        timesteps = self.get_timesteps()  
        if timesteps is not None:
            for t in timesteps:
                outInfo.Append(executive.TIME_STEPS(), t)
            outInfo.Append(executive.TIME_RANGE(), timesteps[0])
            outInfo.Append(executive.TIME_RANGE(), timesteps[-1])
            
            #allocate number of cells that will beadded to the polyline
            #for each timestep we add a cell
            self.pdo.Allocate(self.timeStepsCount ,1)
                    
        return 1
        
 
        
    #RequestData is called for the initial data loading and also when advancing the timesteps
    def RequestData(self, request, inInfoVec, outInfoVec):
  
        data_time = self.get_update_time(outInfoVec.GetInformationObject(0)) 
        
        #the array index of the currently requested time
        currentAnimationTime = np.where(self.time == data_time)
        
        self.get_nc_data(data_time)
        self.drawPolyLine(currentAnimationTime[0][0])
                    
        output = vtk.vtkPolyData.GetData(outInfoVec, 0)
        output.ShallowCopy(self.pdo)
        
        #populate the output with the point and cell data
        for i in self.nc_point_variables:
            output.GetPointData().AddArray(i)
            
        for j in self.nc_cell_variables:
            output.GetCellData().AddArray(j)
 
        if data_time is not None:
            self.pdo.GetInformation().Set(self.pdo.DATA_TIME_STEP(), data_time)
 
        return 1

