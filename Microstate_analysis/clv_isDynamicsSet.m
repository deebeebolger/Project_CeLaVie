function hasDyn = clv_isDynamicsSet(in)
hasDyn = false;
% check if set includes msinfo
if ~isfield(in,'msinfo')
    return;
end
% check if set has FitPar
if ~isfield(in.msinfo, 'FitPar')
    return;
end
% check if FitPar contains Rectify/Normalize parameters
if ~isfield(in.msinfo.FitPar, 'Rectify')
    return;
else
    hasDyn = true;
end
end