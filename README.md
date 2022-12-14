# Estimation of vertex-wise sulcal width maps
Method for the estimation of vertex-wise sulcal width of the brain. 

To execute the method, first you need to compile the travelDepth mex. For this, in Linux console navigate to the mex file and execute:

```
mkdir build
cd build
cmake ..
make
```

Then in MATLAB you can execute the sulcal width estimation algorithm by calling the function EstimateSulcalWidth. This function takes as first parameter the route to the pial file, and as second argument the output where the map will be written (in FreeSurfer format). Optional parameters are allowed such as 

```
'DepthMap' -> Specify a precomputed depth map file
'OutputDepth' -> Specify the output for the depth map file
'DepthStep' -> Parameter to regulate the depth distance between isolines
'DepthThreshold' -> Distance at which the first isoline will be placed
'MaxWidth' -> Maximum width distance allowed to avoid spurious connections
```
## External Dependencies

The repository needs some external dependencies in the MATLAB path to work properly. First of all, you need to clone https://github.com/yasseraleman/Surface_Tools

Then, you need to clone and compile the mex files in https://github.com/alecjacobson/gptoolbox

## TODO

- Add external dependencies as git submodules
- Properly document and comment the code
- Some refactors in the code (split functions, etc.) might be needed
- Create a Docker image to easily run the method
- (FUTURE) Migrate the code from MATLAB to C++/Python
