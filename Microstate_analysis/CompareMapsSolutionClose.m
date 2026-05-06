function CompareMapsSolutionClose(~, ~,fh)
    if isfield(fh.UserData,'CorrMatFig')
        if isvalid(fh.UserData.CorrMatFig)
            delete(fh.UserData.CorrMatFig);
        end
    end

    if fh.UserData.wait
        uiresume(fh);
    else
        delete(fh);
    end
end