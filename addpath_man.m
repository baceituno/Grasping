%adds all of the quadruped code to the path
display('adding code to path');
gurobi_update
root = fileparts(mfilename('fullpath'));
addpath(genpath(root))
