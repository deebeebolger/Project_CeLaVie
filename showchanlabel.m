function output_txt= showchanlabel(~, event, chanlabels)

    idx = event.Position;
    output_txt = {['Participant: ', num2str(idx(1))]; ['Channel index: ', num2str(idx(2))]; ['Channel label: ' chanlabels{idx(2)}]};

end