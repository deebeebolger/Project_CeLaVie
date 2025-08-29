function Filename = Compare_DisplayMDS_MSSolutions(SelectedEEG, nClasses, Filename)

    CompFigHandle = figure('Units', 'normalized', 'Position', [0.2 0.1 0.6 0.8], ...
        'Name', 'Compare microstate maps', 'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'none');

    CompFigHandle.UserData.SelectedEEG = SelectedEEG;
    CompFigHandle.UserData.nClasses = nClasses;
    CompFigHandle.UserData.Filename = [];
    
    CompFigHandle.UserData.CompAxis  = subplot('Position',[0.05 0.13 0.67 0.85],'Parent',CompFigHandle);

    if nClasses == 0
        MinClasses = SelectedEEG.msinfo.ClustPar.MinClasses;
        MaxClasses = SelectedEEG.msinfo.ClustPar.MaxClasses;
        choice = sprintf('%i Classes|', MinClasses:MaxClasses);
        choice(end) = [];
        idx = 1:MaxClasses - MinClasses + 1;
    else
        setIDs = string(1:numel(SelectedEEG));
        CompFigHandle.UserData.setIDs = setIDs;
        CompFigHandle.UserData.setnames = {SelectedEEG.setname};

        choice = strcat(setIDs, ': ', {SelectedEEG.setname});
        idx = 1:numel(SelectedEEG);
    end
   
    obj.Value = idx;
    
    CompareMapsSolutionChanged(obj, [], CompFigHandle);
    
    if isempty(Filename)
        CompFigHandle.UserData.wait = 1;
        viewBtnPos = [0.76 0.09 0.22 0.06];
        uicontrol('style', 'pushbutton', 'string', 'Export shared variances', 'Callback', {@ExportCorrs, CompFigHandle}, 'Units', 'normalized', 'Position', [.76 .02 .22 .06]);
    else
        CompFigHandle.UserData.wait = 0;
        viewBtnPos = [.76 .05 .22 .06];
    end
    uicontrol('style', 'listbox', 'string', choice, 'Value', idx,'min',1,'max',10, ...
        'Callback',{@CompareMapsSolutionChanged,CompFigHandle},'Units','normalized','Position',[0.76 0.17 0.22 0.39]);
    uicontrol('style', 'pushbutton', 'string', 'View shared variances', 'Callback',{@CompareMapsSolutionCorrsToggle,CompFigHandle},'Units','normalized','Position',viewBtnPos);
    uicontrol('style', 'pushbutton' , 'string', '-'   ,  'Callback',{@CompareMapsSolutionScale,CompFigHandle,1.2  },'Units','normalized','Position',[0.225  0.02 0.1 0.05]);
    uicontrol('style', 'pushbutton' , 'string', '+'   ,  'Callback',{@CompareMapsSolutionScale,CompFigHandle,1/1.2},'Units','normalized','Position',[0.335  0.02 0.1 0.05]);
    uicontrol('style', 'pushbutton' , 'string', '>|<' ,  'Callback',{@CompareMapsSolutionScale,CompFigHandle,0},'Units','normalized','Position'    ,[0.445  0.02 0.1 0.05]);
    CompFigHandle.CloseRequestFcn = {@CompareMapsSolutionClose,CompFigHandle};

    if isempty(Filename)
        uiwait(CompFigHandle);
    end

    Filename = CompFigHandle.UserData.Filename;

    % if isvalid(CompFigHandle)
    %     delete(CompFigHandle);
    % end
end