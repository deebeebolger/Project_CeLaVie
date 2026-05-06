
function UpdateCorrTable(fh)
    UserData = fh.UserData;
    if ~isfield(UserData,'CorrMatFig')
        return;
    end
    
    if ~isvalid(UserData.CorrMatFig)
        return;
    end
    UserData.CorrMatFig.Children.Data = UserData.CorrelationTable;
    removeStyle(UserData.CorrMatFig.Children);    
    [row,col] = find(UserData.HitMatrix | eye(size(UserData.HitMatrix,1)));
    tblStyle = uistyle('BackgroundColor',[0.6 1.0 0.6],'FontWeight','bold');
    addStyle(UserData.CorrMatFig.Children,tblStyle,'cell',[row,col]);
end