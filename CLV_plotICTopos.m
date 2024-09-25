function CLV_plotICTopos(EEGin, icrej, ICThresh)

% Subfunction to plot IC topographies isolated by ICLabel.
% It plots 2 figures:
% - A first that shows all ICs with the labels attributed by ICLabel function.
% - A second that presents only those above-threshold ICs marked for
% rejection.
%*****************************************************************************

icClass = string(EEGin.etc.ic_classification.ICLabel.classes);
icClassifed = EEGin.etc.ic_classification.ICLabel.classifications;
icIndx = find(ICThresh(:,1));   % Indices of above-threshold ICs.

assignin('base', 'icrej', icrej);
assignin('base', 'ICThresh', ICThresh)

pop_topoplot(EEGin, 0, EEGin.icachansind, 'IC topographies', 0,'iclabel', 'on')

figIC = figure;
set(figIC, 'Color', [1 1 1], 'Name', 'Above-threshold ICs')
rows = 3;
cols = ceil(length(icrej)/rows);
for counter = 1:length(icrej)

    subplot(rows, cols, counter);
    topoplot(EEGin.icawinv(:, icrej(counter)), EEGin.chanlocs, 'electrodes', 'on');
    title(string(icrej(counter)));

end

end  % end of plotICTopos