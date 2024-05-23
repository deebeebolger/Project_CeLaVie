function plotICTopos(EEGin, icrej, ICThresh)
% Subfunction to plot IC topographies isolated by ICLabel.

icClass = string(EEGin.etc.ic_classification.ICLabel.classes);
icIndx = find(ICThresh(:,1));

pop_topoplot(EEGin, 0)
figure;
for counter = 1:length(icrej)
    subplot(1, length(icrej), counter);
    topoplot(EEGin.icawinv(:, icrej(counter)), EEGin.chanlocs, 'electrodes', 'on');
    title([icClass{icIndx}, '-related Components', string(icrej(counter))]);
end

end  % end of plotICTopos