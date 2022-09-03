# Matlab code: polyvector flow

> Note: final implementation was done in C++, this version is for fast experiments only.

Main file: `polyvector_flow.m`

## minFunc
We use `minFunc` by Mark Schmidt [www.cs.ubc.ca/~schmidtm/Software/minFunc.html](https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html)
 
 1. Download and unzip the content of `minFunc_2012.zip` here: 
  ```
  unzip minFunc_2012.zip
  mv minFunc_2012/* .
  ```
  This will add folders `autoDif/`, `logisticExample/` and `minFunc/`, and files `rosenbrock.m`, `mexAll.m`, `example_minFunc.m`, `example_derivativeCheck.m`.

 2. In Matlab, run this in the current directory: (from `minFunc` instructions)
 ```
>> addpath(genpath(pwd)) % Add all sub-directories to the path
>> mexAll                % Compile mex files (not necessary on all systems)
 ```
