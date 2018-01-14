%adds all of the quadruped code to the path
display('adding code to path');
addpath_drake
root = fileparts(mfilename('fullpath'));
addpath(genpath(root))
