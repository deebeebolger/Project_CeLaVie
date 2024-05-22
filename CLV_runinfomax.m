function [Outeeg, wICA, A, W, IC] = CLV_runinfomax(EEG, RELAX_cfg)

    fprintf('Carry out Independent Components Analysis (ICA) by applying the informax algorithm');
    
    Outeeg = EEG;
    % Maybe run this for the electrodes of type EEG.
    iseeg = ismember({Outeeg.chanlocs.type}, 'EEG');
    ncomps = iseeg;   % Define the number ICs
    [weights, sphere] = runica(Outeeg.data(iseeg,:), 'ncomps', length(ncomps));
    
    %% Copy the IC weigths and sphere information to EEG dataset.
    Outeeg.icaweights = weights;
    Outeeg.icasphere  = sphere;
    Outeeg.icachansind = iseeg;
    Outeeg.icaact = (Outeeg.icaweights*Outeeg.icasphere)*Outeeg.data(iseeg,:); % Calculate the ICA activations.
    Outeeg.icawinv = inv(weights*sphere);

    %% update weight and inverse matrices etc...
    % -----------------------------------------
    if ~isempty(fig), try, close(fig); catch, end; end
    if isempty(Outeeg.icaweights)
        Outeeg.icaweights = pinv(Outeeg.icawinv);
    end
    if isempty(Outeeg.icasphere)
        Outeeg.icasphere  = eye(size(Outeeg.icaweights,2));
    end
    if isempty(Outeeg.icawinv)
        Outeeg.icawinv    = pinv(Outeeg.icaweights*Outeeg.icasphere); % a priori same result as inv
    end

    %% Reorder components by variance
    %  ------------------------------
    meanvar = sum(Outeeg.icawinv.^2).*sum(transpose((Outeeg.icaweights *  Outeeg.icasphere)*Outeeg.data(Outeeg.icachansind,:)).^2)/((length(Outeeg.icachansind)*Outeeg.pnts)-1);
    [~, windex] = sort(meanvar);
    windex = windex(end:-1:1); % order large to small
    meanvar = meanvar(windex);
    Outeeg.icaweights = Outeeg.icaweights(windex,:);
    Outeeg.icawinv    = pinv(Outeeg.icaweights *  Outeeg.icasphere );
    if ~isempty(Outeeg.icaact)
        Outeeg.icaact = Outeeg.icaact(windex,:,:);
    end

    
    
end 