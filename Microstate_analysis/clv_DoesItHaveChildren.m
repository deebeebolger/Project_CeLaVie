function Answer = clv_DoesItHaveChildren(in)
Answer = false;
if ~isfield(in,'msinfo')
    return;
end

if ~isfield(in.msinfo,'children')
    return
else
    Answer = true;
end
end