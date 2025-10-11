function pathSec =  setOSpathSep

if ismac
    pathSec =  '/';
elseif ispc
    pathSec =  '\';
else
    disp('Platform not supported')
end