% adds all of the manipulation code to the path and gurobi (for MACs)
display('adding code to path');
gurobi_update
root = fileparts(mfilename('fullpath'));
addpath(genpath(root))
display('please add gurobi to path')