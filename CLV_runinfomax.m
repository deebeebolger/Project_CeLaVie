function [Outeeg, wIC_all, A, W, IC] = CLV_runinfomax(EEG, RELAX_cfg, fname_curr)

    fprintf('Carry out Independent Components Analysis (ICA) by applying the informax algorithm');
    %% Calculate the rank of the data
    Outeeg = EEG;
    rankeeg = getrank(Outeeg);
    fprintf('The rank of the current dataset is %f.\n', rankeeg);

    %% Maybe ensure that this is run for the electrodes of type EEG only.
    iseeg = ismember({Outeeg.chanlocs.type}, 'EEG');
    ncomps = iseeg;   % Define the number ICs
    [weights, sphere] = runica(Outeeg.data(iseeg,:), 'ncomps', length(ncomps));
    
    %% Copy the IC weigths and sphere information to EEG dataset.
    Outeeg.icaweights = weights;
    Outeeg.icasphere  = sphere;
    Outeeg.icachansind = iseeg;
    W = weights * sphere;
    A = inv(W);
    Outeeg.icawinv = A;
    Outeeg.icaact = (Outeeg.icaweights*Outeeg.icasphere)*Outeeg.data(iseeg,:); % Calculate the ICA activations.

    %% update weight and inverse matrices etc...
    % -----------------------------------------
   
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

    %% Save the data with calculated ICs prior to IC rejection.
    
    fname_ica1  = [fname_curr, '-ICA'];
    Outeeg.setname = fname_ica1;
    saveeeg_ica = fullfile(RELAX_cfg.OutputPath,fname_ica1);
    Outeeg = pop_saveset( Outeeg, saveeeg_ica); 

    %% Identify artifactual ICs with ICLabel.
    %  USE ICLabel to identify components corresponding to artefacts. Then
    %  perform wICA on those components identified.

    Outeeg = iclabel(Outeeg); 

    % IC classes: Brain, Muscle, Eye, Heart, LineNoise, ChannelNoise, Other
    icThreshold    = [0 0.2;0.8 1; 0.8 1; 0.8 1; 0.8 1; 0.8 1; 0 0];
    
    Outeeg = pop_icflag(Outeeg, icThreshold);
    ic2Rej = find(Outeeg.reject.gcompreject);        % Find component/s to reject.
    CLV_plotICTopos(Outeeg, ic2Rej, icThreshold);

    %% Wavelet thresholding of the ICs detected by ICLabel. 
    %  Reference: Neil Baily, 2020. 
    %  Need to define L: level for the stationary wavelet transform.
    %  Default value for L is 5.
    %  Here the data is padded for proper wavelet transform.
    %  The data must be divisible by 2^L. 
    
    if ~isempty(ic2Rej)
        fprintf("Carrying out wavelet thresholding of ICs.\n,..." + ...
            "Begin by padding the data to ensure proper wavelet transform.\n")
       
        L = 5;
        mult = 1;
        wavename = 'coif5';
        IC = Outeeg.icaact; 
        modulus = mod(size(IC,2),2^L); %2^level (level for wavelet)
        if modulus ~=0
            extra = zeros(1,(2^L)-modulus);
        else
            extra = [];
        end

    %% Perform wavelet thresholding.
        fprintf('Perform wavelet thresholding...\n')
        for icnt = 1:numel(ic2Rej)
            if ~isempty(extra)
                sig = [IC(icnt,:),extra]; % pad with zeros
            else
                sig = IC(icnt,:);
            end
            [thresh,sorh,~] = ddencmp('den','wv',sig); % Get automatic threshold value
            thresh = thresh*mult;                      % Multiply threshold by scalar
            swc = swt(sig,L,wavename);                 % Use stationary wavelet transform (SWT) to wavelet transform the ICs
            Y = wthresh(swc,sorh,thresh);              % Threshold the wavelet to remove small values
            wIC(icnt,:) = iswt(Y,wavename);            % Perform inverse wavelet transform to reconstruct a wavelet IC (wIC)
        end
        %% Pad non-artifacted components in the same manner as artifacted components.

        icIndx   = 1:size(IC,1);
        nonartIC = ~ismember(icIndx, ic2Rej); % Find indices of non-artifacted ICs.
        wic_count = 1;
        for icnt2 = 1:numel(nonartIC)
            if nonartIC(icnt2)==1
               wIC_all(icnt2,:) = zeros(1,size(Outeeg.data,2));
               wIC_all(icnt2,:) = [wIC_all(icnt2,:),extra]; % pad with zeros
            elseif nonartIC(icnt2)==0
                wIC_all(icnt2,:) = wIC(wic_count,:);
            end
        end

        %% Plot the wIC and ICs for comparison

        % Remove extra padding
        if ~isempty(extra)
            wIC_all = wIC_all(:,1:end-numel(extra));
        end
        Fs = Outeeg.srate;
        disp('Plotting');
        subplot(3,1,1);
            multisignalplot(IC(ic2Rej,:),Fs,'r');
            title('ICs');
        subplot(3,1,2);
            multisignalplot(wIC_all(ic2Rej,:),Fs,'r');
            title('wICs')
        subplot(3,1,3);
            multisignalplot(IC(ic2Rej,:)-wIC_all(ic2Rej,:),Fs,'r');
            title('Difference (IC - wIC)');

    %% Remove wICA artifact and reconstruct signals. 

    artifacts = A*wIC_all;
    
    % Reshape EEG signal from EEGlab format to channelsxsamples format
    EEG2D = reshape(Outeeg.data, size(Outeeg.data,1), []);
    
    % Subtract out wavelet artifact signal from EEG signal
    wavcleanEEG = EEG2D - artifacts;
    Outeeg.data = wavcleanEEG;


    else
        fprintf(['No artifacted ICs were detected by ICLabel.\n, ...' ...
            'You may need to adjust the icthreshold values.\n']);

    end % end of if isempty(ic2Rej) statement.

    %% ------- Subfunction to calculation the rank of the data -------------
   
    function tmprank2 = getrank(tmpdata)
    
        tmprank = rank(tmpdata);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Here: alternate computation of the rank by Sven Hoffman
        %tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
        covarianceMatrix = cov(tmpdata', 1);
        [E, D] = eig (covarianceMatrix);
        rankTolerance = 1e-7;
        tmprank2=sum (diag (D) > rankTolerance);
        if tmprank ~= tmprank2
            fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because running under Linux 64-bit Matlab\n', tmprank, tmprank2);
            tmprank2 = max(tmprank, tmprank2);
        end

    end

    
    
end 