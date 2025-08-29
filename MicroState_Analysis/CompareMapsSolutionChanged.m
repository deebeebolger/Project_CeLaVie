function CompareMapsSolutionChanged(obj, event, CompFig)
    
    FigAxes = CompFig.UserData;

    MapCollection    = [];
    ColorCollection  = [];
    LabelCollection  = [];
    nClassCollection = [];
    setNumCollection = [];
    CLabelCollection = [];

    if FigAxes.nClasses == 0
        SelectedSolutions = FigAxes.SelectedEEG.msinfo.ClustPar.MinClasses + obj.Value -1;
        for i = SelectedSolutions
            MapCollection     = [MapCollection; FigAxes.SelectedEEG.msinfo.MSMaps(i).Maps];
            ColorCollection   = [ColorCollection; FigAxes.SelectedEEG.msinfo.MSMaps(i).ColorMap];
            for j = 1:i
                LabelCollection   = [LabelCollection FigAxes.SelectedEEG.msinfo.MSMaps(i).Labels(j)];
                nClassCollection  = [nClassCollection i];
                CLabelCollection  = [CLabelCollection,sprintf("%s (%i)",FigAxes.SelectedEEG.msinfo.MSMaps(i).Labels{j},i)];
            end
        end
    else
        SelectedSets = obj.Value;
        for i=SelectedSets
            MapCollection = [MapCollection; FigAxes.SelectedEEG(i).msinfo.MSMaps(FigAxes.nClasses).Maps];
            ColorCollection = [ColorCollection; FigAxes.SelectedEEG(i).msinfo.MSMaps(FigAxes.nClasses).ColorMap];
            for j=1:FigAxes.nClasses
                setNumCollection = [setNumCollection i];
                LabelCollection = [LabelCollection FigAxes.SelectedEEG(i).msinfo.MSMaps(FigAxes.nClasses).Labels(j)];
                CLabelCollection = [CLabelCollection, sprintf("%s (%s)", FigAxes.SelectedEEG(i).msinfo.MSMaps(FigAxes.nClasses).Labels{j}, FigAxes.setnames{i})];
            end
        end
    end    
    MapCollection=double(MapCollection);
    
    ops = statset('Display','off','TolFun',0.01);

    CorrMat = MyCorr(MapCollection').^2;

    pro = mdscale(1-CorrMat,2,'Options',ops,'Criterion','strain');    
  
    cla(FigAxes.CompAxis);
    hold(FigAxes.CompAxis,"on");

    ClassesToDisplay = unique(LabelCollection);
    LegendSubSet = [];

    HitMatrix = false(numel(LabelCollection));
    
%    figure;
%    d = squareform(pdist(pro));
%    plot(d(:),CorrMat(:),'.k')
    
    for i = 1:numel(ClassesToDisplay)
        ItemsToPlot = find(strcmp(ClassesToDisplay{i},LabelCollection));
        nItemsToPlot = numel(ItemsToPlot);
        PlotColor = ColorCollection(ItemsToPlot(1),:);
        
%        clear res;
%        res(:,1) = LabelCollection(ItemsToPlot);
%        res(:,2) = num2cell(nClassCollection(ItemsToPlot))
        
        if nItemsToPlot < 3        
            plot(FigAxes.CompAxis,pro(ItemsToPlot,1),pro(ItemsToPlot,2),'-k');
        else
            xy = pro(ItemsToPlot,:);
            Hull = convhull(xy(:,1),xy(:,2));
            patch(FigAxes.CompAxis,xy(Hull,1),xy(Hull,2),PlotColor,'FaceAlpha',.2);
        end
        ph = plot(FigAxes.CompAxis,pro(ItemsToPlot,1),pro(ItemsToPlot,2),'ok','MarkerFaceColor',PlotColor,'MarkerSize',8);        
        for j = 1:nItemsToPlot
            idx = ItemsToPlot(j);
            if FigAxes.nClasses == 0                
                txt = sprintf(" %i",nClassCollection(idx));                
            else
                txt = [ ' ' FigAxes.setIDs(setNumCollection(idx))];
            end
            text(FigAxes.CompAxis,pro(idx,1),pro(idx,2),txt,'HorizontalAlignment','center', ...
                    'VerticalAlignment','Bottom','Interpreter','none','FontSize', 11, 'FontWeight', 'bold');
            HitMatrix(idx,ItemsToPlot) = true;
        end
        LegendSubSet = [LegendSubSet,ph];
    end
    
    legend(LegendSubSet,ClassesToDisplay,'Interpreter','none','Position',[0.76,0.58,0.22,0.39]);

    axis(FigAxes.CompAxis,'equal');
    axis(FigAxes.CompAxis,'tight');    
    lim = max(abs([FigAxes.CompAxis.XLim FigAxes.CompAxis.YLim]));
    lim = lim * 1.2;
    FigAxes.OrgLim = lim;

    FigAxes.CompAxis.XLim = [-lim lim];
    FigAxes.CompAxis.YLim = [-lim lim];
    
    FigAxes.CorrelationTable = array2table(CorrMat * 100,'VariableNames',CLabelCollection,'RowNames',CLabelCollection);
    FigAxes.HitMatrix = HitMatrix;

    UpdateCorrTable(CompFig);

    CompFig.UserData = FigAxes;    
end