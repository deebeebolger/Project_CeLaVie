
function T = CLV_LoadMetaData(xlspath, xlsfname)


xls_fullpath = fullfile(xlspath, xlsfname);
sheet = 'Posttest';
T = readtable(xls_fullpath, sheet='Posttest', ReadVariableNames=true);

end 