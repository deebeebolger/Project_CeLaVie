function isPublished = clv_isPublishedSet(in, templateNames)
isPublished = false;
if isempty(in.setname)
    return;
end

if matches(in.setname, templateNames)
    isPublished = true;
end
end