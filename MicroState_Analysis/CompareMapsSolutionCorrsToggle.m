function CompareMapsSolutionCorrsToggle(~, ~,fh)
    
    if isfield(fh.UserData,'CorrMatFig')
        if isvalid(fh.UserData.CorrMatFig)
            delete(fh.UserData.CorrMatFig);
            return;
        end
    end

    fh.UserData.CorrMatFig = uifigure;

    fh.UserData.CorrMatFig.WindowStyle = "normal";
    fh.UserData.CorrMatFig.Name = "Percent shared variance";
    uitable('Unit','normalized','Position',[0.02 0.02 0.96 0.96],'Parent',fh.UserData.CorrMatFig);
    UpdateCorrTable(fh);
end