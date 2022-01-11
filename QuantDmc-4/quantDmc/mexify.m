%compiles the C functions to MEX

currentPath = pwd;

aFile = 'quantBiDmcMulti.m';
newPath = which(aFile);
if length(newPath) == 0
    error(sprintf('Could not find %s in the search path\nAdd QuantDmc to the search path',upper(aFile)));
end
newPath = strrep(newPath,aFile,'');
cd(newPath)

fprintf('Compiling....')

try
    mex -v clib/quantBiDmc.c clib/mmatrix.c
    fprintf('success\n');
    cd(currentPath)
catch
    fprintf('Error, failed to compile\n');
    cd(currentPath)
end