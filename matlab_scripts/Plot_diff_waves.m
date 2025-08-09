% Script_check_diff_waves

    % Uses the output from Script9_GrandAverage_ERPs.m
    % Uses the output from Permutations_and_analysis.m

    % Plots the calculated ERP difference waveforms to check their compatibility with the
    % results from the ERP CORE publication
    % Calculates Cohen's dz effect size by averaging amplitude difference
    % from 0 uV using individual difference waveforms

    % External dependencies:
    addpath(genpath('./eeglab'));
    addpath(genpath('./erplab'));

close all; clearvars;
 
SUB = {'1','2','3','4','5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40'};
Nsub = length(SUB);
COMP = {'N170', 'MMN', 'N2pc', 'N400','P3', 'LRP','ERN'};

% Specify component specific parameters; list needs to match with order of COMP above!
spec_channels = {'26', '20','9','13','14','5','20'};
diff_bins = {'5','3','1','3','3','1','7'};
measure_windows = {[110 150], [125 225], [200 275], [300 500], [300 600], [-100 0], [0 100]};


% Directory with ERP Core grand average files of all components
grandavDIR = fullfile('Preprocessed_ERP_Core','All_components_files', 'ERP_Core_grand_averages');
% make figure directory if it doesn't exist yet
figDIR = fullfile('Preprocessed_ERP_Core','All_components_files', 'matlab_diff_wave_figures');
    if ~exist(figDIR, 'dir')
        mkdir(figDIR);
    end

excludeMap = containers.Map();
excludeMap('N170') = {'1','5','16'};
excludeMap('MMN') = {'7'};
excludeMap('N2pc') = {'7','9','10','12','28'};
excludeMap('P3') = {'6','9','10','30','35','40'};
excludeMap('N400') = {'40'};
excludeMap('LRP') = {'6','30','40'};
excludeMap('ERN') = {'5','6','30','40'};

% Initialize table to save effect sizes (mean amplitude)
effect_size_res = array2table(nan(length(SUB), length(COMP)), 'VariableNames', COMP);

%% Make a "data-check" plot of diff waves obtained and those from ERP CORE
for c = 1:length(COMP)
   
    currentDIR = fullfile('Preprocessed_ERP_Core','All_components_files', COMP{c});
    dataDIR = fullfile(currentDIR, 'data');
    resultsDIR = fullfile(currentDIR, 'perm_results'); 

    % Initialize table to save group-level waveform and CIs (all components)
    mat = load(fullfile(resultsDIR,'2_results.mat'),'res');
    Xf      = mat.res.Xf; % vector of time points (ms)
    group_vars = {'Xf','c1_av','c1_ci_low', 'c1_ci_upp', 'c2_av', 'c2_ci_low', 'c2_ci_upp', 'diff_av', 'diff_ci_low','diff_ci_upp'};
    group_av_res = array2table(nan(length(Xf), 10), 'VariableNames', group_vars);
    group_av_res.Xf = Xf';

    % Initialize matrices to temporarily store wave forms from all participants in
    all_c1_spec = nan(length(SUB), length(Xf));
    all_c2_spec = nan(length(SUB), length(Xf));
    all_diff_spec = nan(length(SUB), length(Xf));


% get index of participants to exclude
  exclude = excludeMap(COMP{c});
  final_Nsub = Nsub - length(exclude);


% get index of bins and electrode of interest
channelIdx = str2double(spec_channels{c});
binIdx = str2double(diff_bins{c});
measurement_window = measure_windows{c};

% Get component specific baseline window and epoch window for x-axis
    if strcmp(COMP{c}, 'LRP') % baseline from -800 to -600 ms 
        baseWin = [-800 -600];
        xscale = [-800.0 200.0];
    elseif strcmp(COMP{c}, 'ERN') % baseline from -400 to -200 ms
        baseWin = [-400 -200];
        xscale = [-600.0 400.0];
    else
        baseWin = [-199 0]; % prestimulus phase is bsaeline for all other components
        xscale = [-200.0 800.0];
    end

% make a plot for each component of all individual diff waveforms
figure("Visible","off"); hold on;
h = gobjects(0); % empty array for handles
legend_labels = {};

for i = 1:Nsub 
    % skip excluded participants
    if ismember(SUB{i}, exclude)
            disp(['Skipping subject ' SUB{i} ' for component ' COMP{c}]);
            continue;
    end

    % get ERP diff waves calculated from Permutations_and_analysis.m
    mat = load(fullfile(resultsDIR,[SUB{i} '_results.mat']),'res');
    diff_spec      = mat.res.diff_spec; 
    av_c1_spec     = mat.res.av_c1_spec;
    av_c2_spec     = mat.res.av_c2_spec;

    % save to all-subjects matrix
    all_c1_spec(i, :) = av_c1_spec; 
    all_c2_spec(i, :) = av_c2_spec; 
    all_diff_spec(i, :) = diff_spec; 

    
    % Baseline corrected difference wave
    baseline = Xf >= baseWin(1) & Xf <= baseWin(2);
    diff_spec_bn = diff_spec - mean(diff_spec(baseline));

    % Calculate effect size  of the difference in mean amplitude from 0 μV
    measure = Xf >= measurement_window(1) & Xf <= measurement_window(2);
    effect_size_res{i, c} = mean(diff_spec_bn(measure),'omitna');
   
     
    % plot individual diff waveforms of all subjects
    h(end+1) = plot(Xf, diff_spec);
    legend_labels{end+1} = ['participant ' SUB{i}];

end % subject loop


% calculate average and confidence intervals for group-level averages of condition 1
group_av_res.c1_av = mean(all_c1_spec, 1,'omitnan')';
group_sem_c1 = std(all_c1_spec, 0, 1, 'omitnan') / sqrt(sum(~isnan(all_c1_spec),1));
group_av_res.c1_ci_upp = group_av_res.c1_av + 1.96 * group_sem_c1;
group_av_res.c1_ci_low = group_av_res.c1_av - 1.96 * group_sem_c1;

% calculate average and confidence intervals for group-level averages of condition 2
group_av_res.c2_av = mean(all_c2_spec, 1,'omitnan')';
group_sem_c2 = std(all_c2_spec, 0, 1, 'omitnan') / sqrt(sum(~isnan(all_c2_spec),1));
group_av_res.c2_ci_upp = group_av_res.c2_av + 1.96 * group_sem_c2;
group_av_res.c2_ci_low = group_av_res.c2_av - 1.96 * group_sem_c2;

% calculate group level average of difference waveform
group_av_res.diff_av = mean(all_diff_spec, 1,'omitnan')';
group_sem_diff = std(all_diff_spec, 0, 1, 'omitnan') / sqrt(sum(~isnan(all_diff_spec),1));
group_av_res.diff_ci_upp = group_av_res.diff_av + 1.96 * group_sem_diff;
group_av_res.diff_ci_low = group_av_res.diff_av - 1.96 * group_sem_diff;

    % finish figure with individual difference waves
    yline(0, 'k', 'LineWidth', 1); % black zero reference line
    plot(Xf, group_av_res.diff_av, 'k', 'LineWidth', 2); % grand average line
    title(['Individual Difference Waveforms for ' COMP{c}]);
    xlabel('Time (ms)');
    ylabel('Amplitude (µV)');
    xlim(xscale);
    saveas(gcf, fullfile(figDIR, [COMP{c} '_indiv_diff.png']));
    hold off;

% GRAND AVERAGE DIFF WAVEFORM and ERP Core grand average

% get ERP grand average diff waves from ERP CORE Script Script9_GrandAverage_ERPs.m
ERP = pop_loaderp('filename', ['GA_' COMP{c} '_erp_ar_diff_waves.erp'], 'filepath', grandavDIR);  
ERP_Core_diff_erp = squeeze(ERP.bindata(channelIdx, :, binIdx));

% Plot average difference wave
figure('Visible','off'); hold on;
plot(Xf, group_av_res.diff_av, 'k', 'LineWidth', 2); 
yline(0, 'k', 'LineWidth', 1); % black zero reference line
% Plot SEM shading
fill_x = [Xf, fliplr(Xf)];  % x for shading polygon
fill_y = [group_av_res.diff_ci_low', fliplr(group_av_res.diff_ci_upp')];
fill(fill_x, fill_y, 'k', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
% Plot ERP CORE difference wave
plot(Xf, ERP_Core_diff_erp, 'b', 'LineWidth', 1);
xlabel('Time (ms)');
ylabel('Amplitude (µV)');
xlim(xscale);
title(['Group Average Difference Waveform of ' COMP{c}]);
legend({'Difference Wave', 'Confidence Interval', 'ERP CORE data'}, 'Location', 'Best');
saveas(gcf, fullfile(figDIR, [COMP{c} '_group_diff.png']));
hold off;

% Save table with group average waveforms and CIs
writetable(group_av_res, fullfile(currentDIR, [COMP{c} '_group_av_CIs.xlsx']));

end % component loop

% Save table with effect size calculations
writetable(effect_size_res, fullfile('Preprocessed_ERP_Core', 'Effect_size_results.xlsx'));


