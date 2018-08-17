% adds all of the manipulation code with Drake and Gurobi to the path
display('adding code to path');
addpath_drake
gurobi_update
root = fileparts(mfilename('fullpath'));
addpath(genpath(root))