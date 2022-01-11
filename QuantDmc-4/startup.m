addpath([pwd '/quantDmc']);

if exist('quantBiDmc') ~= 3
    warning(sprintf('QUANTBIDMC MEX file is not available for this platform\nRun MEXIFY to compile.\n'));
    mexify
end