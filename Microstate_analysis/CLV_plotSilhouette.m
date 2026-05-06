function meanSilh = CLV_plotSilhouette(MapsUsed, mIndx, rows, cols, count, nClust)
    
   
    subplot(rows, cols, count)
    [silh1, h] = silhouette(MapsUsed,mIndx,'correlation');
    xlabel('Silhouette Value')
    ylabel('Cluster')

    meanSilh = mean(silh1);
    titleStr = nClust + "-Cluster Solution ( " + meanSilh + " )";
    title(titleStr)

end

