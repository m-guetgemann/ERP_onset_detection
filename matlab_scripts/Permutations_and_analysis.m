% Script of Meta Gütgemann for analysis of ERP Core data 

% Input: operates on preprocessed ERP Core data (output of Preprocessing_ERP_Core.m)

% Output: obtains univariate t^2 and multivariate HT^2 from empirical (H1) and permutated (H0) data
    % saves all variables within the res structure as perm_results.mat files 
    % in the perm_results folder seperately for each ERP Component
    % also saves low-pass filtered ERP parent & difference waveforms

% External Dependencies
    % Uses the limo toolbox https://github.com/LIMO-EEG-Toolbox/limo_tools
    % Uses highdim toolbox https://github.com/brian-lau/highdim/blob/master/%2Bdiff/hotell2.m 
addpath(genpath('./functions'));
addpath(genpath('./eeglab'));
addpath(genpath('./erplab'));

close all; clearvars;

% List of ERP Components to process
COMP = {'N170', 'MMN','N2pc', 'N400','P3', 'LRP','ERN'};

% ERP Core selected electrodes of interest and bins (Script12_Measure_ERPs.m)
diff_bins = {'5','3','1','3','3','1','7'};
selected_channels = {'26', '20','14','5','20'}; % matches order of COMP  
    % N170; PO8 (chan 26)
    % MMN; FCz (chan 20)
    % N2pc; PO7/PO8 (chan 9/26)
    % P3; Pz (chan 13)
    % N400; CPz (chan 14)
    % LRP; C3/C4 (chan 5/22)
    % ERN; FCz (chan 20)

% Full subjects list
SUB = {'1', '2', '3', '4','5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40'};
Nsub = length(SUB);

% Specify Parameters
ath = 0.05; % set alpha level
pth = 1 - ath;  % obtain permutation threshold
Nperm = 2000; % set number of permutations, ERP Core uses 2500 for clustering
rng(1); 

% Define list of participants to exclude per component
    % this list follows the processing steps of the ERP Core publication
excludeMap = containers.Map();
excludeMap('N170') = {'1','5','16'};
excludeMap('MMN') = {'7'};
excludeMap('N2pc') = {'7','9','10','12','28'};
excludeMap('P3') = {'6','9','10','30','35','40'};
excludeMap('N400') = {'40'};
excludeMap('LRP') = {'6','30','40'};
excludeMap('ERN') = {'5','6','30','40'};

% open EEG lab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
scalpIdx = [1:28];

%**********************************************************************************************************************************************************************

%% Loop through each ERP Component listed in COMP
for c = 1:length(COMP)

    % Specify directories for the current component
    currentDIR = fullfile('Preprocessed_ERP_Core','All_components_files', COMP{c});
    dataDIR = fullfile(currentDIR, 'data');
    

    % Create 'figures', 'perm_results' and 'seeds' directories
    figDIR = fullfile(currentDIR, 'figures');
    if ~exist(figDIR, 'dir')
        mkdir(figDIR);
    end
    resultsDIR = fullfile(currentDIR, 'perm_results');
    if ~exist(resultsDIR, 'dir')
        mkdir(resultsDIR);
    end
    seedsDIR = fullfile(currentDIR, 'seeds');
    if ~exist(seedsDIR, 'dir')
        mkdir(seedsDIR);
    end

    % Get component specific parameters
    channelIdx = str2double(selected_channels{c});
    diff_binIdx = str2double(diff_bins{c});

% Get excluded subjects for this component
    if isKey(excludeMap, COMP{c})
        exclude = excludeMap(COMP{c});
    else
        exclude = {};
    end

% Loop through each subject listed in SUB
for i = 1:Nsub
    
    if ismember(SUB{i}, exclude)
            disp(['Skipping subject ' SUB{i} ' for component ' COMP{c}]);
            continue;
        end
    
disp(['Now processing component ' COMP{c} ' of participant ' SUB{i}]);


% Load the preprocessed data file outputted from Preprocessing_ERP_Core.m
EEG = pop_loadset( 'filename', [SUB{i} '_prepro.set'], 'filepath', dataDIR);
Xf = EEG.times;
% Get component specific electrde of interest label
elecLabel = EEG.chanlocs(channelIdx).labels;
% Get component specific baseline window
    if strcmp(COMP{c}, 'LRP') % baseline from -800 to -600 ms
        baseline = Xf >= -800 & Xf <= -600;
        baseWin = [-800 -600];
    elseif strcmp(COMP{c}, 'ERN') % baseline from -400 to -200 ms
        baseline = Xf >= -400 & Xf <= -200;
        baseWin = [-400 -200];
    else
        baseline = Xf < 0; % baseline is prestimulus period for all other components
        baseWin = [-199 0];
    end
    
% reject epochs flagged with artifacts
rejEpochs = EEG.reject.rejmanual; 
EEG = pop_rejepoch(EEG, rejEpochs, 0); 

% Low-pass filter and baseline normalize EEG data
    EEG = pop_basicfilter(EEG,[], 'Filter', 'lowpass', 'Design', 'butter','Cutoff', 20,'Order', 8);            

% Exclude rereferencing channels 
EEG = pop_select(EEG, 'channel', scalpIdx); 

% Seperate data from appropriate bin labels (EEG.event.bini)
    % LRP and N2pc contrasts are obtained from contra- minus ipsi-lateral sites
    %LH = [1,2,3,4,5,6,7,8,9,10,11];      
    %RH = [15,17,18,19,22,23,24,25,26,27,28];
switch COMP{c}
    % select main electrode of interest (as specified by ERP Core)
    case 'N2pc'
        LH = 9; RH = 26; % PO7/PO8
        % Extract the first bin number from each epoch's eventbini 
        bini_all = arrayfun(@(x) x.eventbini{1}(1), EEG.epoch);
    case 'LRP'
        LH = 5; RH = 22; % C3/C4
        % Extract relevant bin from each epoch (LRP saves 2 bini per epoch)
        bini_all = arrayfun(@(e) EEG.event(e.event(min(end,2))).bini(1), EEG.epoch);
end

if strcmp(COMP{c}, 'N2pc') || strcmp(COMP{c}, 'LRP')
    % bin 1; left target (N2pc) or left response (LRP)
    % bin 2; right target (N2pc) or right response (LRP)
    b1_RH = EEG.data(RH,:,bini_all == 1);  
    b2_LH = EEG.data(LH,:,bini_all == 2);  
    b1_LH = EEG.data(LH,:,bini_all == 1);  
    b2_RH = EEG.data(RH,:,bini_all == 2);  
end
switch COMP{c}
    case 'N2pc'
    % N2pc;  contra = (b1@RH + b2@LH)/2 - ipsi = (b1@LH + b2@RH)/2
    cond1_data = cat(3, b1_RH, b2_LH); % [chan x time x trials] 
    cond2_data   = cat(3, b1_LH, b2_RH);
    case 'LRP'
    % LRP; ipsi = contra = (b1@RH + b2@LH)/2 - (b1@LH + b2@RH)/2 
    cond1_data = cat(3, b1_RH, b2_LH); % [chan x time x trials] 
    cond2_data   = cat(3, b1_LH, b2_RH);
    otherwise
    % Default condition: non-lateralized ERP contrast at all electrodes
    cond1EEG = pop_selectevent(EEG, 'bini', 1);  % experimental condition
    cond1_data = cond1EEG.data;
    cond2EEG = pop_selectevent(EEG, 'bini', 2); % control condition
    cond2_data = cond2EEG.data;
end

% Get trial number
[Ne, Nf, Ntc1] = size(cond1_data);
[~, ~, Ntc2] = size(cond2_data);

% Trial-average waveform by condition
av_c1 = mean(cond1_data, 3);   
av_c2 = mean(cond2_data, 3);
av_c1_bn = av_c1 - hd(av_c1(baseline)); 
av_c2_bn = av_c2 - hd(av_c1(baseline));
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    av_c1_spec = av_c1_bn(channelIdx,:);  
    av_c2_spec = av_c2_bn(channelIdx,:);
else
    av_c1_spec = av_c1;
    av_c2_spec = av_c2;
end

% Difference waveform 
diff_erp = av_c1 - av_c2;     % LRP and N2pc; contra- ipsilateral
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    diff_spec = diff_erp(channelIdx, :);
else
    diff_spec = diff_erp; % N2pc and LRP variables only have one channel anyway
end

% SEM across trials
sem_c1av = std(cond1_data, 0, 3) / sqrt(Ntc1);
sem_c2av = std(cond2_data, 0, 3) / sqrt(Ntc2);

% 95% confidence intervals
t_val_c1 = tinv(1 - ath/2, Ntc1 - 1);
t_val_c2 = tinv(1 - ath/2, Ntc2 - 1);
ci_lower_c1av = av_c1 - t_val_c1 .* sem_c1av;
ci_upper_c1av = av_c1 + t_val_c1 .* sem_c1av;
ci_lower_c2av = av_c2 - t_val_c2 .* sem_c2av;
ci_upper_c2av = av_c2 + t_val_c2 .* sem_c2av;

% Save standard deviation parent waveforms by condition
std_c1 = std(cond1_data,[], 3); 
std_c2 = std(cond2_data,[], 3);

% Confidence interval for the difference waveform, for independent samples
    % report trial-level variability! 
sem_diff = sqrt(sem_c1av.^2 + sem_c2av.^2);
Ndf = min(Ntc1, Ntc2) - 1; % degrees of freedom chosen from smaller trial set
t_val_diff = tinv(1 - ath/2, Ndf);
ci_lower_diff = diff_erp - t_val_diff .* sem_diff;
ci_upper_diff = diff_erp + t_val_diff .* sem_diff;


% Empirical t-test using limo toolbox

[m, dfe, ci, sd, n, tval, pval, tcrit] = limo_ttest(2,cond1_data,cond2_data,ath);
                            % output tval has dimensions electrodes x timepoints
t2 = tval.^2;
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    [maxt2, maxIdx] = max(t2, [], 1); % selects maximum electrode at each timepoint
    max_tval = tval(sub2ind(size(tval), maxIdx, 1:length(maxIdx))); % use Idx for each 
    max_tcrit = tcrit(sub2ind(size(tcrit), maxIdx, 1:length(maxIdx)));
    max_pval = pval(sub2ind(size(pval), maxIdx, 1:length(maxIdx)));

    spec_m = m(channelIdx,:);
    spec_dfe = dfe(channelIdx,:);
    spec_ci = ci(channelIdx,:);
    spec_tval = tval(channelIdx,:);
    spec_t2 = t2(channelIdx,:);
    spec_tcrit = tcrit(channelIdx,:);
    spec_pval = pval(channelIdx,:);
else
    spec_m = m;
    spec_dfe = dfe;
    spec_ci = ci;
    spec_tval = tval;
    spec_t2 = t2;
    spec_pval = pval;
    spec_tcrit = tcrit;
end

% plot t^2 
figname = ['participant ',SUB{i}];
figure('Name', figname,'NumberTitle','off','Color','white', 'Visible','off');
xline(0, 'k', 'LineWidth', 1); hold on;  % vertical black line at time zero
h1 = plot(Xf, t2, 'Color','k', 'LineWidth', 0.5);
hl1 = h1(1);
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    h2 = plot(Xf, maxt2, 'r', 'LineWidth', 1);   % maxt2 virtual electrode, red
    h3 = plot(Xf, t2(channelIdx,:), 'Color','g','LineWidth',1); % ERP Core electrode of interest
end
set(gca, 'Layer', 'top');
xlabel('Time in ms');
ylabel('t^2');
title(['Participant ' SUB{i} ' : univariate t^2 at all electrodes']);
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    legend([hl1 h2 h3], {'t^2 (all electrodes)', 'max t^2', ['electrode ' elecLabel]});
    else
    legend(hl1, ['t^2 from electrode ' elecLabel]);
end
ax = gca;
ax.FontSize = 14; 
ylim([min(0, min(ylim)), max(ylim)]);
filename = fullfile(figDIR, [SUB{i} '_fig_t2.png']);
saveas(gcf, filename);

% Empirical MANOVA using hotell2 function from highdim toolbox
HT2_allNe = zeros(1, Nf);       % preallocate for all electrodes (one T2 per frequency/time)
%HT2_subsetNe = zeros(1, Nf);    % preallocate for subset of electrodes
HT2_gradient = zeros(Ne, Nf);   % preallocate for gradients (channels × freq/time)

% using all electrodes
 for F = 1:Nf 
    x = squeeze(cond1_data(:, F, :))'; % now [trials x channels]
    y = squeeze(cond2_data(:, F, :))'; 
    [pval, T2] = hotell2(x,y);
    HT2_allNe(F) = T2; % extract Hotelling's T2
 end

 % using ERP gradients as features for MANOVA
        lin_cond1 = zeros(size(cond1_data)); % (format Ne x Nf x Nt)
        lin_cond2 = zeros(size(cond2_data));
        for E = 1:Ne
            for T = 1:Ntc1
                lin_cond1(E, :, T) = gradient(squeeze(cond1_data(E,:,T)),2);
            end
            for T = 1:Ntc2
                lin_cond2(E, :, T) = gradient(squeeze(cond2_data(E, :,T)),2);
            end
        end
  % compute MANOVA with gradient as features
  for E = 1:Ne 
    for F = 1:Nf
        temp_cond1 = squeeze(cond1_data(E, F, :));    % Ntc1 trials 
        temp_cond2 = squeeze(cond2_data(E, F, :));  % Ntc2 trials
        temp_lin_cond1 = squeeze(lin_cond1(E,F,:));
        temp_lin_cond2 = squeeze(lin_cond2(E,F,:));
        [pval, T2] = hotell2([temp_cond1, temp_lin_cond1],[temp_cond2, temp_lin_cond2]);
        HT2_gradient(E,F) = T2; % Ne x Nf
    end
  end
  if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    maxHT2_gradient = max(HT2_gradient, [], 1);
    spec_HT2_gradient = HT2_gradient(channelIdx,:);
    else
      spec_HT2_gradient = HT2_gradient(:,:);
  end

% Plot Hotelling's T^2 from MANOVA Results
figname = ['participant ', SUB{i}];
figure('Name', figname,'NumberTitle','off','Color','white', 'Visible','off');
xline(0, 'k', 'LineWidth', 1); hold on;  % vertical black line at time zero
h1 = plot(Xf, HT2_gradient', 'Color','k', 'LineWidth', 0.5); 
hl1 = h1(1);
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    h2 = plot(Xf, maxHT2_gradient, 'Color','r', 'LineWidth', 1);
    h3 = plot(Xf, HT2_gradient(channelIdx,:), 'Color', 'g', 'LineWidth', 1); % ERP Core electrode of interest
end
xlabel('Time in ms');
ylabel('Hotellings T^2');
title(['Participant ' SUB{i} ' : multivariate HT^2 on amplitude and gradient']);
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    legend([hl1 h2 h3], {'HT^2 (all electrodes)', 'max HT^2', ['electrode ' elecLabel]});
    else
    legend(hl1, ['HT^2 from electrode ' elecLabel]);
end
ax = gca;
ax.FontSize = 14; 
ylim([min(0, min(ylim)), max(ylim)]);
filename = fullfile(figDIR, [SUB{i} '_fig_HT2.png']);
saveas(gcf, filename);

% Obtain Permutation estimates 

% Create Permutation index 
indx_perm = zeros(Nperm, Ntc1 + Ntc2);
for perm_iter = 1:Nperm % Nperm permutations in 100 rows
    indx_perm(perm_iter,:) = randperm(Ntc1 + Ntc2); 
end

% Concatenate data from Condition 1 and 2
allEEG = cat(3,cond1_data, cond2_data); % along 3rd dimension of trials

% Initialize arrays for permutation estimates
tval_perm = zeros(Ne,Nf,Nperm); % channels x timepoints x permutations
t2_perm = zeros(Ne, Nf, Nperm);
spec_t2_perm = zeros(Nf,Nperm);
%subsetHT2_perm = zeros(Nf, Nperm);        % for MANOVA subset electrodes
gradientHT2_perm = zeros(Ne, Nf, Nperm);  % for MANOVA gradients
spec_HT2_perm = zeros(Nf, Nperm);  
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    maxt2_perm = zeros(Nf,Nperm);
    maxHT2_perm = zeros(Nf, Nperm);           % max gradient HT2 across electrodes
end

% initialize arrays for baseline normalized values
t2_perm_bn = zeros(Ne, Nf, Nperm);
spec_t2_perm_bn = zeros(Nf,Nperm);
%subsetHT2_perm_bn = zeros(Nf, Nperm);
gradientHT2_perm_bn = zeros(Ne, Nf, Nperm);
spec_HT2_perm_bn = zeros(Nf, Nperm); 
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    maxt2_perm_bn = zeros(Nf, Nperm);
    maxHT2_perm_bn = zeros(Nf, Nperm);
end

for perm_iter = 1:Nperm
    perm_trials = indx_perm(perm_iter,:);

    perm_c1 = allEEG(:,:,perm_trials(1:Ntc1));
    perm_c2 = allEEG(:,:,perm_trials(Ntc1+1:Ntc1+Ntc2));

    if perm_iter == 1
    disp(['Now computing P' SUB{i} ' permutation no. ' num2str(perm_iter) ' of ' COMP{c}])
    end

    if rem(perm_iter, 100) == 0
    disp(['Now computing P' SUB{i} ' permutation no. ' num2str(perm_iter) ' of ' COMP{c}])
    end

    % t-test
    [~, ~, ~, ~, ~, tval, ~] = limo_ttest(2,perm_c1,perm_c2,ath);
    t2_temp = tval.^2;
    tval_perm(:,:,perm_iter) = tval;
    t2_perm(:,:,perm_iter) = t2_temp;
    
    if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
        spec_t2_perm(:,perm_iter) = t2_temp(channelIdx, :);
        maxt2_perm(:,perm_iter) = max(squeeze(t2_temp), [], 1);
        else
        spec_t2_perm(:,perm_iter) = t2_temp(:, :);
    end

    % baseline normalized t2
    t2_perm_bn(:,:,perm_iter) = t2_temp - hd(t2_temp(baseline));
    if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
        spec_t2_perm_bn(:,perm_iter) = t2_perm_bn(channelIdx, :, perm_iter);
        maxt2_perm_bn(:, perm_iter) = max(squeeze(t2_perm_bn(:,:,perm_iter)), [], 1);
        else
        spec_t2_perm_bn(:,perm_iter) = t2_perm_bn(:, :, perm_iter);
    end
    
    % manova at electrodes independently, using gradients as features, to get H(0) max HT2
        lin_perm_c1 = zeros(size(perm_c1)); % (format Ne x Nf x Nt)
        lin_perm_c2 = zeros(size(perm_c2));
        for E = 1:Ne
            for T = 1:Ntc1
                lin_perm_c1(E, :, T) = gradient(squeeze(perm_c1(E,:,T)),2);
            end % end Ntc1
            for T = 1:Ntc2
                lin_perm_c2(E, :, T) = gradient(squeeze(perm_c2(E, :,T)),2);
            end % end Ntc2
        end % end E

  % compute MANOVA permutation data (gradient as features)
  for E = 1:Ne 
    for F = 1:Nf
        temp_perm_c1 = squeeze(perm_c1(E, F, :));    % Ntc1 trials 
        temp_perm_c2 = squeeze(perm_c2(E, F, :));  % Ntc2 trials
        temp_lin_perm_c1 = squeeze(lin_perm_c1(E,F,:));
        temp_lin_perm_c2 = squeeze(lin_perm_c2(E,F,:));
        [pval, T2] = hotell2([temp_perm_c1, temp_lin_perm_c1],[temp_perm_c2, temp_lin_perm_c2]);
        gradientHT2_perm(E,F,perm_iter) = T2; % Ne x Nf x Nperm
    end % end Nf
  end % end Ne
  if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
     spec_HT2_perm(:, perm_iter) = squeeze(gradientHT2_perm(channelIdx, :, perm_iter));
     maxHT2_perm(:, perm_iter) = max(squeeze(gradientHT2_perm(:,:, perm_iter)), [], 1);
  else
      spec_HT2_perm(:, perm_iter) = squeeze(gradientHT2_perm(:, :, perm_iter));
  end

  % baseline normalized HT2 gradient
  baseline_tmp = hd(reshape(gradientHT2_perm(:, baseline, perm_iter), [], 1));
  gradientHT2_perm_bn(:,:,perm_iter) = gradientHT2_perm(:,:,perm_iter) - baseline_tmp;
  if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    spec_HT2_perm_bn(:, perm_iter) = squeeze(gradientHT2_perm_bn(channelIdx, :, perm_iter));
    maxHT2_perm_bn(:, perm_iter) = max(squeeze(gradientHT2_perm_bn(:,:,perm_iter)), [], 1);
  else
    spec_HT2_perm_bn(:, perm_iter) = squeeze(gradientHT2_perm_bn(:, :, perm_iter));
  end

end % permutations loop

% Plots to check data from first permutation (without baseline normalization)
figure('Name', 'tval_perm','Visible','off'); plot(Xf,squeeze(tval_perm(:, :, 1)));

%figure('Name', 'subsetHT2_perm','Visible','off'); plot(Xf, squeeze(subsetHT2_perm(:,1)));
figure('Name', 'gradientHT2_perm','Visible','off'); plot(Xf, squeeze(gradientHT2_perm(:, :, 1)));


% Plot t-test and MANOVA results (baseline normalized results)

% Baseline normalisation -- subtract 50th quantile of baseline T2 and HT2
    % for permutated values already implemented in perm_iter loop above
    t2_bn = t2 - hd(t2(baseline));
    %HT2_subsetNe_bn = HT2_subsetNe - hd(HT2_subsetNe(baseline));
    HT2_gradient_bn = HT2_gradient - hd(HT2_gradient(baseline));
    if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
        spec_t2_bn = t2_bn(channelIdx, :);
        maxt2_bn = max(t2_bn, [],1);
        maxHT2_gradient_bn = max(HT2_gradient_bn,[], 1);
        spec_HT2_bn = HT2_gradient_bn(channelIdx, :);
    else
        spec_t2_bn = t2_bn(:, :);
        spec_HT2_bn = HT2_gradient_bn(:, :);
    end

% Plot permutation estimates of max t2 and empirical max t2
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    figname = ['Participant ',SUB{i}];
    figure('Name', figname,'NumberTitle','off','Color','white', 'Visible','off');
    xline(0, 'k', 'LineWidth', 1); hold on;  % vertical black line at time zero
    h1 = plot(Xf, maxt2_perm_bn, 'Color', [.7 .7 .7], 'LineWidth', 0.5); % grey; 1000 permutations of maxt2
    hl1 = h1(1);
    h2 = plot(Xf, maxt2_bn, 'r', 'LineWidth', 1); hold on; % red; empirical maxt2
    set(gca, 'Layer', 'top');
    xlabel('Time in ms');
    ylabel('t^2');
    title(['Participant ' SUB{i} ' : empirical and permutated max t^2']);
    legend([hl1 h2], {'permutated max t^2', 'empirical max t^2'});
    ax = gca;
    ax.FontSize = 14; 
    ylim([min(0, min(ylim)), max(ylim)]);
    filename = fullfile(figDIR, [SUB{i} '_fig_maxt2_perm.png']);
    saveas(gcf, filename);
end

% Plot permutation and empirical t2 statistic from electrode of interest
figname = ['Participant ',SUB{i}];
figure('Name', figname,'NumberTitle','off','Color','white', 'Visible','off');
xline(0, 'k', 'LineWidth', 1); hold on;  % vertical black line at time zero
h1 = plot(Xf, spec_t2_perm_bn, 'Color', [.7 .7 .7], 'LineWidth', 0.5); % grey; 1000 permutations of maxt2
hl1 = h1(1);
h2 = plot(Xf, spec_t2_bn, 'r', 'LineWidth', 1); hold on; % red; empirical maxt2
set(gca, 'Layer', 'top');
xlabel('Time in ms');
ylabel('t^2');
title(['Participant ' SUB{i} ' : t^2 from electrode of interest']);
legend([hl1 h2], {'permutated t^2', 'empirical t^2'});
ax = gca;
ax.FontSize = 14; 
ylim([min(0, min(ylim)), max(ylim)]);
filename = fullfile(figDIR, [SUB{i} '_fig_spec_t2_perm.png']);
saveas(gcf, filename);

% Plot permutation estimates of max HT2 virtual electrode and empirical max HT2
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
figname = ['Participant ', SUB{i}];
figure('Name', figname,'NumberTitle','off','Color','white', 'Visible','off');
xline(0, 'Color', 'k', 'LineWidth', 1); hold on;  % vertical black line at time zero
h1 = plot(Xf, maxHT2_perm_bn, 'Color', [.7 .7 .7], 'LineWidth', 0.5); % grey; 1000 permutations of maxHT2
hl1 = h1(1);
h2 = plot(Xf, maxHT2_gradient_bn, 'r', 'LineWidth', 1); hold on; % red; empirical maxHT2
set(gca, 'Layer', 'top');
xlabel('Time in ms');
ylabel('Hotellings T^2');
title(['Participant ' SUB{i} ' : empirical and permutated max HT^2']);
legend([hl1 h2], {'permutated max HT^2', 'empirical max HT^2'});
ax = gca;
ax.FontSize = 14; 
ylim([min(0, min(ylim)), max(ylim)]);
filename = fullfile(figDIR, [SUB{i} '_fig_maxHT2_perm.png']);
saveas(gcf, filename);
end

% Plot permutation and empirical HT2 statistic from electrode of interest
figname = ['Participant ', SUB{i}];
figure('Name', figname,'NumberTitle','off','Color','white', 'Visible','off');
xline(0, 'Color', 'k', 'LineWidth', 1); hold on;  % vertical black line at time zero
h1 = plot(Xf, spec_HT2_perm_bn, 'Color', [.7 .7 .7], 'LineWidth', 0.5); % grey; 1000 permutations of maxHT2
hl1 = h1(1);
h2 = plot(Xf, spec_HT2_bn, 'r', 'LineWidth', 1); hold on; % red; empirical maxHT2
set(gca, 'Layer', 'top');
xlabel('Time in ms');
ylabel('Hotellings T^2');
title(['Participant ' SUB{i} ' : HT^2 from electrode of interest']);
legend([hl1 h2], {'permutated HT^2', 'empirical HT^2'});
ax = gca;
ax.FontSize = 14; 
ylim([min(0, min(ylim)), max(ylim)]);
filename = fullfile(figDIR, [SUB{i} '_fig_spec_HT2_perm.png']);
saveas(gcf, filename);

% Save results 
    % Variables are commented out to save computation space and time
    % but can be un-commented to save for different analysis purposes

% save seed for participant and component (in case script crashes)
seed = rng;
save(fullfile(seedsDIR, [SUB{i},'_', COMP{c}, '_perm_seed.mat']), 'seed');

% save t2 and HT2 results from original and permutated data
res = [];
res.Xf = Xf;

% save erp parent waveforms
res.av_c1 = av_c1; % at all electrodes here!
res.av_c2 = av_c2; 
res.sem_c1av = sem_c1av;
res.sem_c2av = sem_c2av;
%res.std_c1 = std_c1;
%res.std_c2 = std_c2;

% save these variables only for the specified electrode of interest
res.av_c1_spec = av_c1_spec;
res.av_c2_spec = av_c2_spec;
res.ci_lower_c1av = ci_lower_c1av; 
res.ci_upper_c1av = ci_upper_c1av;
res.ci_lower_c2av = ci_lower_c2av;
res.ci_upper_c2av = ci_upper_c2av;

% save erp difference waveforms
res.diff_erp = diff_erp; % at all electrodes here!
res.sem_diff = sem_diff;
res.diff_spec = diff_spec; % only electrode of interest
res.ci_lower_diff = ci_lower_diff;
res.ci_upper_diff = ci_upper_diff;

% Save univariate t-test results from electrode of interest
res.spec_m = spec_m;
res.spec_dfe = spec_dfe;
res.spec_ci = spec_ci;
res.spec_tval = spec_tval;
res.spec_t2 = spec_t2;
res.spec_pval = spec_pval;
res.spec_tcrit = spec_tcrit;

% Save baseline-normalized t2 and HT2 time-series
res.spec_t2_bn = spec_t2_bn;      
res.spec_HT2_bn = spec_HT2_bn;
%res.HT2_gradient = HT2_gradient;       % slows down script, not needed
%res.spec_HT2_gradient = spec_HT2_gradient; 
%res.t2_bn = t2_bn;                    
%res.HT2_gradient_bn = HT2_gradient_bn; 

% Save baseline-normalized permutation results
res.spec_t2_perm_bn = spec_t2_perm_bn;
res.spec_HT2_perm_bn = spec_HT2_perm_bn;
%res.spec_t2_perm = spec_t2_perm;           % slows down script, not needed
%res.spec_HT2_perm = spec_HT2_perm;
%res.t2_perm_bn = t2_perm_bn;            
%res.gradientHT2_perm_bn = gradientHT2_perm_bn; 

% Save max t2 and max HT2 variables only for given components
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    res.maxt2 = maxt2;
    res.max_tval = max_tval;
    res.max_tcrit = max_tcrit;
    res.max_pval = max_pval;

    res.maxt2_bn = maxt2_bn;   
    res.maxHT2_gradient_bn = maxHT2_gradient_bn; 
    res.maxt2_perm_bn = maxt2_perm_bn;
    res.maxHT2_perm_bn = maxHT2_perm_bn;

    %res.maxHT2_gradient = maxHT2_gradient;   % not needed, slows down script
    %res.maxt2_perm = maxt2_perm;
    %res.maxHT2_perm = maxHT2_perm;
end

save(fullfile(resultsDIR, [SUB{i}, '_results.mat']), 'res');
disp('The seed and permutations estimates have been saved as a .mat file');

end % participant loop

disp(['Finished analysing all subjects of component: ' COMP{c}]);

end % component loop

disp('Finished analysing all components');



