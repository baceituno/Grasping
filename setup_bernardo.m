% adds all of the manipulation code with Drake and Gurobi to the path
display('adding code to path');
% addpath_drake
current_fldr = pwd;
cd /Library/gurobi810/mac64/examples/matlab
cd ../..
cd matlab
gurobi_setup
cd (current_fldr)
root = fileparts(mfilename('fullpath'));
addpath(genpath(root))