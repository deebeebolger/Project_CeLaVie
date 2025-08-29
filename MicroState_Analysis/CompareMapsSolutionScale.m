function CompareMapsSolutionScale(~, ~,fh, Scaling)
    UserData = fh.UserData;
    FigAxes = UserData;
    if Scaling == 0
        lim = UserData.OrgLim;
        FigAxes.CompAxis.XLim = [-lim lim];
        FigAxes.CompAxis.YLim = [-lim lim];

    else
        FigAxes.CompAxis.XLim = FigAxes.CompAxis.XLim * Scaling;
        FigAxes.CompAxis.YLim = FigAxes.CompAxis.YLim * Scaling;
    end
end