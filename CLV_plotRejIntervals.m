function CLV_plotRejIntervals()
%% Programmer: Deirdre BOLGER Date: 08-10-2024
%  Simple function permitting the user to visualise those time intervals that were automatically rejected during the detection of extreme data.
%  It is necessary to load in two datasets, one after the other (not ideal,
%  I will fix this):
%  - first dataset is the dataset just after resampling and before bad data rejection.
%  - second dataset is the final dataset after interpolation.
%  The time intervals that were rejected will be highlighted on a
%  multichannel plot of the resampled data (that has not yet been cleaned).
%*********************************************************************************

rechoix = 1;
prompt = {'Load another dataset (Yes/No):'};

fnameAll = {};
fpathAll = {};

while rechoix == 1

    [filename, filepath] = uigetfile({'*.bdf; *.set' 'BDF or eeglab set file';'*.bdf' 'BDF'; '*.set' 'eeglab set file'}, 'Choose a BDF or .set data file ', 'multiselect', 'on');
    [parentdir, ~,~] = fileparts(filepath);
    [~,name,ext] = fileparts(filename);
    main_dir = fileparts(parentdir);
    fnameAll = [fnameAll; {filename}];
    fpathAll = [fpathAll; {filepath}];
    
    dlgtitle = 'Input';
    fieldsize = [1 45];
    definput = {'Yes'};
    answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

    if strcmp(answer, 'Yes')
        rechoix = 1;
    else
        rechoix = 0;
    end
    
end 

dataIn = cellfun(@(x,y) pop_loadset(x,y), fnameAll, fpathAll, 'UniformOutput',false); 
datarej = dataIn{2,1}.RELAX.ExtremelyBadPeriodsForDeletion;

winrej = [datarej, repmat([1 1 0],size(datarej,1), 1), ones(size(datarej,1),128)];
eegplot(dataIn{1,1}.data, 'winlength',20, 'winrej', winrej);

end
