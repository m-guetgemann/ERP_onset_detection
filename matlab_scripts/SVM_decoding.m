% Decoding of single-participant ERP data using SVM
    % using multivariate data across all electrodes
    
% Input: preprocessed output from ERP CORE Script #6: Artifact_Rejection.m
% Output: SVM decoding accuracy time-series results
    % follows SVM code by Carrasco et al. (2024)
    % P2b_Best_Prep.m 
    % P2b_Decoding.m

% External Dependencies:
addpath(genpath('./eeglab'));
addpath(genpath('./erplab'));

close all; clearvars;

% Specify ERP components
COMP = {'N170', 'MMN', 'N2pc', 'N400', 'LRP','ERN'};

% Specify subjects 
SUB = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40'};	
Nsub = length(SUB);

% Specify SVM parameters
bin_labels = [1 2]; 
prep_filename = 'prep';
results_filename = 'SVM';
channels = [1:28];
iter = 100;
crossfolds = 4; % Ntrials/nCrossblocks = Ntrials(ERP)


% exclude participants from analysis
excludeMap = containers.Map();
excludeMap('N170') = {'1','5','16'};
excludeMap('MMN') = {'7'};
excludeMap('N2pc') = {'7','9','10','12','28'};
excludeMap('P3') = {'6','9','10','30','35','40'};
excludeMap('N400') = {'40'};
excludeMap('LRP') = {'6','30','40'};
excludeMap('ERN') = {'5','6','30','40'};

% open EEG lab
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% Loop through each component in COMP
for c = 1:length(COMP)

     % Specify directories for the current component
    currentDIR = fullfile('Preprocessed_ERP_Core','All_components_files', COMP{c});
    dataDIR = fullfile(currentDIR, 'data');
  
% Location of folder to save the prepped .best files
prepDIR = fullfile(currentDIR,'prep_for_SVM');
if ~exist(prepDIR,'dir')
    mkdir(prepDIR);
end
% Location of folder to save the decoding results
resultsDIR = fullfile(currentDIR,'SVM_results');
if ~exist(resultsDIR,'dir')
    mkdir(resultsDIR);
end
% Location of folder to save the decoding plots
figDIR = fullfile(currentDIR,'SVM_figures');
if ~exist(figDIR,'dir')
    mkdir(figDIR);
end

% Get component specific epoch window 
    if strcmp(COMP{c}, 'LRP') % epoch from -800 to 200 ms 
       epochWin = [-800 200]; 
    elseif strcmp(COMP{c}, 'ERN') % epoch from -600 to 400 ms 
        epochWin = [-600 400]; 
    else
        epochWin = [-200 800]; 
    end

% Get excluded subjects for this component
    if isKey(excludeMap, COMP{c})
        exclude = excludeMap(COMP{c});
    else
        exclude = {};
    end
%% Loop through each subject listed in SUB

for i = 1:Nsub
    
    if ismember(SUB{i}, exclude)
            disp(['Skipping subject ' SUB{i} ' for component ' COMP{c}]);
            continue;
        end
    
disp(['Now processing component ' COMP{c} ' of participant ' SUB{i}]);
    
    % Load the epoched and artifact rejected EEG data file outputted from Script #6 in .set EEGLAB file format
    EEG = pop_loadset( 'filename', [SUB{i} '_prepro.set'], 'filepath', dataDIR);
    EEG = pop_basicfilter( EEG,  1:35, 'Cutoff',  20, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  8 ); 
        % following Carrasco et al. (2024)

    % Create bin epoched single trial (BEST) structures 
    BEST = pop_extractbest(EEG,'Bins', bin_labels);
    BEST = pop_savemybest(BEST,'filename',[SUB{i} '_' prep_filename], ...
            'filepath',prepDIR, 'overwriteatmenu','on');

    % Run SVM decoding analysis with specified channels and parameters
    MVPC = pop_decoding(BEST, 'Classes', bin_labels, 'Channels', channels, 'nIter', iter, 'nCrossblocks', crossfolds, ...
    'DecodeTimes', epochWin, 'Decode_Every_Npoint', 1, 'EqualizeTrials', 'classes', ...
    'Method', 'SVM', 'ParCompute', 'on');
        % output MVPC.times
        % output MVPC.Accuracy
        % output average_score (average decoding accuracy over iter/folds)
        % output stderror (of decoding accuracy over iter/folds)
        % output raw_preditions (raw classifier output

    % Save SVM decoding results
    MVPC  = pop_savemymvpc(MVPC,'mvpcname',[SUB{i} '_' results_filename], 'filename',[SUB{i} '_' results_filename '.mvpc'], ...
        'filepath', resultsDIR, 'overwriteatmenu', 'on', 'warning', 'off'); % will overwrite existing results!

%  Plot AUC over time

figname = ['Participant ', SUB{i}];
figure('Name', figname,'NumberTitle','off','Color','white','Visible','off');
plot(MVPC.times, MVPC.average_score); hold on;
fill_between = @(x, y1, y2, color) fill([x fliplr(x)], [y1 fliplr(y2)], color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill_between(MVPC.times, MVPC.average_score - MVPC.stderror, MVPC.average_score + MVPC.stderror, 'b'); % shaded std error
xlabel('Time (ms)');
xline(0, 'k', 'LineWidth', 1);
xlim(epochWin);
ylabel('Decoding Accuracy');
title(['Participant ' SUB{i} ' : Decoding Accuracy over Time']);
filename = fullfile(figDIR, [SUB{i} '_SVM_accuracy.png']);
saveas(gcf, filename);

end % subject loop

disp(['Finished analysing all subjects of component: ' COMP{c}]);

end % component loop

disp('Finished analysing all components');
