% Plot results from Permutations_and_analysis.m

addpath(genpath('./eeglab'));
addpath(genpath('./erplab'));

close all; clearvars;

% List of ERP Components to process
COMP = {'N170','MMN','N2pc', 'N400', 'P3', 'LRP', 'ERN'};

% ERP Core selected electrodes of interest and bins (Script12_Measure_ERPs.m)
diff_bins = {'5','3','1', '3', '3','1','7'};
selected_channels = {'26', '20','9', '13','14','5','20'}; % matches order of COMP  
all_elecLabels = {'PO8','FCz','PO7','Pz','CPz','C3','FCz'};
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
Nperm = 2000; % set number of permutations
rng(1); % set seed

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
    % c = 1;
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
    elecLabel = all_elecLabels{c};

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

% i = 1;

% get results calculated from Permutations_and_analysis.m
mat = load(fullfile(resultsDIR,[SUB{i} '_results.mat']),'res');
Xf      = mat.res.Xf; 





if isfield(mat.res, 't2_bn')
spec_t2_perm_bn     = mat.res.t2_perm_bn(channelIdx, :, :);
spec_t2_bn          = mat.res.t2_bn(channelIdx, :);
spec_HT2_perm_bn    = mat.res.gradientHT2_perm_bn(channelIdx, :, :);
spec_HT2_bn         = mat.res.HT2_gradient_bn(channelIdx, :);
else 
spec_t2_bn                  = mat.res.spec_t2_bn;
spec_t2_perm_bn              = mat.res.spec_t2_perm_bn;
spec_t2_perm_bn             = mat.res.spec_HT2_perm_bn;
spec_HT2_perm_bn        = mat.res.spec_HT2_perm_bn;
spec_HT2_bn                 = mat.res.spec_HT2_bn;
end

if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
maxt2_perm_bn          = mat.res.maxt2_perm_bn;
maxt2_perm_bn       = mat.res.maxt2_perm_bn;
maxt2_bn            = mat.res.maxt2_bn;
maxHT2_perm_bn      = mat.res.maxHT2_perm_bn;
maxHT2_gradient_bn  = mat.res.maxHT2_gradient_bn;
end

% plot t^2 
figname = ['participant ',SUB{i}];
figure('Name', figname,'NumberTitle','off','Color','white', 'Visible','off');
xline(0, 'k', 'LineWidth', 1); hold on;  % vertical black line at time zero
h1 = plot(Xf, spec_t2_bn, 'Color','b', 'LineWidth', 1);
hl1 = h1(1);
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    h2 = plot(Xf, maxt2_bn, 'r', 'LineWidth', 1);   % maxt2 virtual electrode, red
end
set(gca, 'Layer', 'top');
xlabel('Time in ms');
ylabel('t^2');
title([COMP{c} ' | Participant ' SUB{i} ' : univariate t^2']);
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    legend([hl1 h2], {['electrode ' elecLabel], 'max t^2'});
    else
    legend([ 't^2 from electrode ' elecLabel]);
end
ax = gca;
ax.FontSize = 14; 
ylim([min(0, min(ylim)), max(ylim)]);
filename = fullfile(figDIR, [SUB{i} '_fig_t2.png']);
saveas(gcf, filename);

% Plot Hotelling's T^2 from MANOVA Results
figname = ['participant ', SUB{i}];
figure('Name', figname,'NumberTitle','off','Color','white', 'Visible','off');
xline(0, 'k', 'LineWidth', 1); hold on;  % vertical black line at time zero
h1 = plot(Xf, spec_HT2_bn', 'Color','b', 'LineWidth', 1); 
hl1 = h1(1);
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    h2 = plot(Xf, maxHT2_gradient_bn, 'Color','r', 'LineWidth', 1);
end
xlabel('Time in ms');
ylabel('Hotellings T^2');
title([COMP{c} ' | Participant ' SUB{i} ' : HT^2 (amplitude & gradient)']);
if ~strcmp(COMP{c}, 'N2pc') && ~strcmp(COMP{c}, 'LRP')
    legend([hl1 h2], {['electrode ' elecLabel],'max HT^2'});
    else
    legend(hl1, ['HT^2 from electrode ' elecLabel]);
end
ax = gca;
ax.FontSize = 14; 
ylim([min(0, min(ylim)), max(ylim)]);
filename = fullfile(figDIR, [SUB{i} '_fig_HT2.png']);
saveas(gcf, filename);


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
    title([COMP{c} ' | Participant ' SUB{i} ' : empirical and permutated max t^2']);
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
h1 = plot(Xf, squeeze(spec_t2_perm_bn), 'Color', [.7 .7 .7], 'LineWidth', 0.5); % grey; 1000 permutations of maxt2
hl1 = h1(1);
h2 = plot(Xf, spec_t2_bn, 'r', 'LineWidth', 1); hold on; % red; empirical maxt2
set(gca, 'Layer', 'top');
xlabel('Time in ms');
ylabel('t^2');
title([COMP{c} ' | Participant ' SUB{i} ' : t^2 from electrode of interest']);
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
title([COMP{c} ' | Participant ' SUB{i} ' : empirical and permutated max HT^2']);
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
h1 = plot(Xf, squeeze(spec_HT2_perm_bn), 'Color', [.7 .7 .7], 'LineWidth', 0.5); % grey; 1000 permutations of maxHT2
hl1 = h1(1);
h2 = plot(Xf, spec_HT2_bn, 'r', 'LineWidth', 1); hold on; % red; empirical maxHT2
set(gca, 'Layer', 'top');
xlabel('Time in ms');
ylabel('Hotellings T^2');
title([COMP{c} ' | Participant ' SUB{i} ' : HT^2 from electrode of interest']);
legend([hl1 h2], {'permutated HT^2', 'empirical HT^2'});
ax = gca;
ax.FontSize = 14; 
ylim([min(0, min(ylim)), max(ylim)]);
filename = fullfile(figDIR, [SUB{i} '_fig_spec_HT2_perm.png']);
saveas(gcf, filename);

end % participant loop

disp(['Finished analysing all subjects of component: ' COMP{c}]);

end % component loop

disp('Finished analysing all components');



