
function hasMS = clv_hasMicrostates(in)
hasMS = false;

% check if set includes msinfo
if ~isfield(in,'msinfo')
    return;
end

% check if msinfo is empty
if isempty(in.msinfo)
    return;
else
    hasMS = true;
end
end