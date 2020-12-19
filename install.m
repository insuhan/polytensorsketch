fprintf("Compiling MEX functions...\n");
cd ./mexfunctions
compile
cd ..

fprintf("Adding paths...\n");
addpath(genpath('.'));
fprintf("Done!\n");