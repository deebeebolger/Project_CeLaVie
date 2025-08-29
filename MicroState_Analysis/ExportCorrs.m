function ExportCorrs(obj, ~, fh)
    CorrTable = fh.UserData.CorrelationTable;
    CorrTable.Properties.RowNames = CorrTable.Properties.VariableNames;

    [FName, PName, idx] = uiputfile({'*.csv', 'Comma separated file'; '*.txt', 'Tab delimited file'; '*.xlsx', 'Excel file'; '*.mat', 'Matlab Table'}, 'Save shared variance matrix');
    if FName == 0
        return;
    end
    
    Filename = fullfile(PName, FName);
    if idx < 4
        writetable(CorrTable, Filename, 'WriteRowNames', true);
    else
        save(Filename, 'CorrTable');
    end

    fh.UserData.Filename = [fh.UserData.Filename {Filename}];
end