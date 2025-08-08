% Script by Meta Gütgemann to obtain preprocessed EEG data 
    % this script is a combination of all #1 - #6 preprocessing scripts
    % from the ERP CORE Project and the credit goes to it's authors: 
    % Kappenman, E., Farrens, J., Zhang, W., Stewart, A. X., & Luck, S. J. (2020). ERP CORE: An Open Resource for Human Event-Related Potential Research [Preprint]. PsyArXiv. https://doi.org/10.31234/osf.io/4azqm

% Runs over all 7 components of the ERP Core 
% Operates on individual subject data (N = 40 participants)
% Uses the output from Script #3: Run_ICA.m (semi-continuous EEG data file containing ICA weights from Script #3 + loads list of ICA component(s) from the ICA_Components_' COMP{c} '.xlsx Excel file)
% (Note; if ICA weights were re-computed on the data, the component(s) to remove will need to be updated in the Excel file to match the new components (see Script #3: Run_ICA.m for further details))

% ROUGH OVERVIEW OF PROCESSING STEPS;
% 1. Removes ICA component(s) from the EEG
% 2. Create an Event List containing a record of all event codes and their timing
% 3. Assign events to bins using Binlister (and event markers specific to each component)
% 4. Epoch the EEG and perform baseline correction (specific to each component)
% 5. Interpolate bad channels listed in Excel file Interpolate_Channels_' COMP{c} '.xls
% 6. Perform artifact rejection to remove noisy segments of EEG segments containing eyeblinks or eye movements
% (uses individual subject's parameters that are listed in the corresponding Excel file for each artifact)

% External dependencies:
addpath(genpath('./eeglab'));
addpath(genpath('./erplab'));

close all; clearvars;

% Specify location of folder containing this script and all raw data folders
    % download ERP CORE raw data folders from https://osf.io/thsqg/
    % files are slow to download on flashdrive; faster on computer harddrive
baseDIR = 'C:\Users\meta\Desktop\Data\Matlab-projects\ERP-Core'


% Specify a new folder to save the final pre-processed data files
%preprocessed_dataDIR = 'E:\Preprocessed_ERP_Core';
preprocessed_dataDIR = fullfile(baseDIR, 'preprocessed_ERN')

% List of ERP Components to process
COMP = {'ERN_raw'};


%**********************************************************************************************************************************************************************

%Loop through each ERP Component listed in COMP
for c = 1:length(COMP)
    currentComp = COMP{c};

%Location of the main study directory for the current component
DIR = [fullfile(baseDIR, currentComp) '_raw'];

%Location of folder that contains associated processing files for each component
Current_File_Path = fullfile(DIR, 'EEG_ERP_Processing');

%List of subjects to process, based on the name of the folder that contains that subject's data
SUB = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40'};    

%Load the Excel file with the list of channels to interpolate for each subject 
[ndata1, text1, alldata1] = xlsread([Current_File_Path filesep 'Interpolate_Channels_' COMP{c}]);

%Load the Excel file with the list of thresholds and parameters for identifying C.R.A.P. with the simple voltage threshold algorithm for each subject 
[ndata2, text2, alldata2] = xlsread([Current_File_Path filesep 'AR_Parameters_for_SVT_CRAP_' COMP{c}]);

%Load the Excel file with the list of thresholds and parameters for identifying C.R.A.P. with the moving window peak-to-peak algorithm for each subject 
[ndata3, text3, alldata3] = xlsread([Current_File_Path filesep 'AR_Parameters_for_MW_CRAP_' COMP{c}]);


%Load the Excel file with the list of thresholds and parameters for identifying any uncorrected horizontal eye movements (using the ICA-corrected HEOG signal) with the step like algorithm for each subject 
[ndata6 text6, alldata6] = xlsread([Current_File_Path filesep 'AR_Parameters_for_SL_HEOG_' COMP{c}]);

if ~strcmp(COMP{c}, 'MMN') && ~strcmp(COMP{c}, 'LRP') && ~strcmp(COMP{c}, 'ERN');
 % not for MMN and LRP and ERN component
%Load the Excel file with the list of thresholds and parameters for identifying eyeblinks during the stimulus presentation period (using the original non-ICA corrected VEOG signal) with the moving window peak-to-peak algorithm for each subject 
[ndata4, text4, alldata4] = xlsread([Current_File_Path filesep 'AR_Parameters_for_MW_Blinks_' COMP{c}]);

%Load the Excel file with the list of thresholds and parameters for identifying horizontal eye movements during the stimulus presentation period (using the original non-ICA corrected HEOG signal) with the step like algorithm for each subject 
[ndata5 text5, alldata5] = xlsread([Current_File_Path filesep 'AR_Parameters_for_SL_HEOG_Stim_Pres_' COMP{c}]);
end

%**********************************************************************************************************************************************************************

%Loop through each subject listed in SUB
for i = 1:length(SUB)

    %Open EEGLAB and ERPLAB Toolboxes  
    [ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
    
    %Define subject path based on study directory and subject ID of current subject
    Subject_Path = [DIR filesep SUB{i} filesep];

    % SCRIPT #4 from ERP Core; Script4_Remove_ICA_Components.m

    %Load the continuous EEG data file containing the ICA weights outputted from Script #3 in .set EEGLAB file format
    if strcmp(COMP{c}, 'MMN')
    EEG = pop_loadset('filename',[SUB{i} '_' COMP{c} '_ds_reref_ucbip_hpfilt_ica_weighted.set'],'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[SUB{i} '_' COMP{c} '_ds_reref_ucbip_hpfilt_ica_weighted'], 'gui','off'); 
    else
    EEG = pop_loadset('filename',[SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_weighted.set'],'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',[SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_weighted'], 'gui','off'); 
    end

    %Load list of ICA component(s) corresponding to ocular artifacts from Excel file ICA_Components_N170.xlsx
    [ndata, text, alldata] = xlsread([Current_File_Path filesep 'ICA_Components_' COMP{c}]); 
    MaxNumComponents = size(alldata, 2);
        for j = 1:length(alldata)
            if isequal(SUB{i}, num2str(alldata{j,1}));
                NumComponents = 0;
                for k = 2:MaxNumComponents
                    if ~isnan(alldata{j,k});
                        NumComponents = NumComponents+1;
                    end
                    Components = [alldata{j,(2:(NumComponents+1))}];
                end
            end
        end

    %Perform ocular correction by removing the ICA component(s) specified above
    EEG = pop_subcomp( EEG, [Components], 0);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',[SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr'],'savenew', [Subject_Path SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr.set'],'gui','off'); 
    
    if strcmp(COMP{c}, 'LRP') && strcmp(COMP{c}, 'ERN')

    %Create a bipolar HEOG channel (HEOG_left minus HEOG_right) and a bipolar VEOG channel (VEOG_lower minus FP2) from the ICA corrected data; the original uncorrected HEOG and VEOG channels are retained for later artifact detection procedures
    EEG = pop_eegchanoperator( EEG, [Current_File_Path filesep 'Add_Corrected_Bipolars_' COMP{c} '.txt']);

    else

    %Create a bipolar HEOG channel (HEOG_left minus HEOG_right) and a bipolar VEOG channel (VEOG_lower minus FP2) from the ICA corrected data; the original uncorrected HEOG and VEOG channels are retained for later artifact detection procedures
    EEG = pop_eegchanoperator( EEG, [Current_File_Path filesep 'Add_Corrected_Bipolars_' COMP{c} '.txt'], 'Saveas', 'off', 'History', 'script');

    end 

    %Add channel location information corresponding to the 3-D coordinates of the electrodes based on 10-10 International System site locations
    EEG = pop_chanedit(EEG, 'lookup',[Current_File_Path filesep 'standard-10-5-cap385.elp']);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip'], 'savenew', [Subject_Path SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip.set'], 'gui', 'off'); 

    % SCRIPT 5 from ERP-Core; Script5_Elist_Bin_Epoch.m
 
    %Load the semi-continuous ICA-corrected EEG data file outputted from Script #4 in .set EEGLAB file format
    EEG = pop_loadset( 'filename', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip.set'], 'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip'], 'gui', 'off'); 

    %Create EEG Event List containing a record of all event codes and their timing
    EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, 'Eventlist', [Subject_Path SUB{i} '_' COMP{c} '_Eventlist.txt'] ); 
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist'], 'savenew', [Subject_Path SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist.set'], 'gui', 'off');

    %Assign events to bins with Binlister; an individual trial may be assigned to more than one bin (bin assignments can be reviewed in each subject's COMP{c} _Eventlist_Bins.txt file)
    EEG  = pop_binlister( EEG , 'BDF', [Current_File_Path filesep 'BDF_' COMP{c} '.txt'], 'ExportEL', [Subject_Path SUB{i} '_' COMP{c} '_Eventlist_Bins.txt'], 'IndexEL',  1, 'SendEL2', 'EEG&Text', 'UpdateEEG', 'on', 'Voutput', 'EEG' );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins'], 'savenew', [Subject_Path SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins.set'], 'gui', 'off'); 

    if strcmp(COMP{c}, 'LRP') 
        %Epoch the EEG into 1-second segments time-locked to the response (from -800 ms to 200 ms) and perform baseline correction using the average activity from -800 ms to -600 ms 
        EEG = pop_epochbin( EEG , [-800.0  200.0],  [-800.0  -600.0]);
    elseif strcmp(COMP{c}, 'ERN')
        %Epoch the EEG into 1-second segments time-locked to the response (from -600 ms to 400 ms) and perform baseline correction using the average activity from -400 ms to -200 ms 
        EEG = pop_epochbin( EEG , [-600.0  400.0],  [-400.0  -200.0]);
    else 
        %Epoch the EEG into 1-second segments time-locked to the response (from -200 ms to 800 ms) and perform baseline correction using the average activity from -200 ms to 0 ms 
        EEG = pop_epochbin( EEG , [-200.0  800.0],  [-200.0  0.0]);
    end 

    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch'], 'savenew', [Subject_Path SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch.set'], 'gui', 'off');
    close all;

    % SCRIPT 6 from ERP-Core; Script6_Artifact_Rejection.m
    %Load the epoched EEG data file outputted from Script #5 in .set EEGLAB file format
    EEG = pop_loadset( 'filename', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch.set'], 'filepath', Subject_Path);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch'], 'gui', 'off'); 

    %Interpolate channel(s) specified in Excel file Interpolate_Channels_' COMP{c} '.xls; any channels without channel locations (e.g., the eye channels) should not be included in the interpolation process and are listed in ignored channels
    %EEG channels that will later be used for measurement of the ERPs should not be interpolated
    ignored_channels = [29 30 31 32 33 34 35];        
    DimensionsOfFile1 = size(alldata1);
    for j = 1:DimensionsOfFile1(1);
        if isequal(SUB{i},num2str(alldata1{j,1}));
           badchans = (alldata1{j,2});
           if ~isequal(badchans,'none') | ~isempty(badchans)
           	  if ~isnumeric(badchans)
                 badchans = str2num(badchans);
              end
              EEG  = pop_erplabInterpolateElectrodes( EEG , 'displayEEG',  0, 'ignoreChannels',  ignored_channels, 'interpolationMethod', 'spherical', 'replaceChannels', badchans);
           end
           [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch_interp'], 'savenew', [Subject_Path SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch_interp.set'], 'gui', 'off'); 
        end
    end

    %Identify segments of EEG with C.R.A.P. artifacts using the simple voltage threshold algorithm with the parameters in the Excel file for this subject
    DimensionsOfFile2 = size(alldata2);
    for j = 1:DimensionsOfFile2(1)
        if isequal(SUB{i},num2str(alldata2{j,1}));
            if isequal(alldata2{j,2}, 'default')
                Channels = 1:31;
            else
                Channels = str2num(alldata2{j,2});
            end
            ThresholdMinimum = alldata2{j,3};
            ThresholdMaximum = alldata2{j,4};
            TimeWindowMinimum = alldata2{j,5};
            TimeWindowMaximum = alldata2{j,6};
        end
    end

    EEG  = pop_artextval( EEG , 'Channel',  Channels, 'Flag', [1 2], 'Threshold', [ThresholdMinimum ThresholdMaximum], 'Twindow', [TimeWindowMinimum  TimeWindowMaximum] ); 
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch_interp_SVT'], 'gui', 'off'); 

    %Identify segments of EEG with C.R.A.P. artifacts using the moving window peak-to-peak algorithm with the parameters in the Excel file for this subject
    DimensionsOfFile3 = size(alldata3);
    for j = 1:DimensionsOfFile3(1)
        if isequal(SUB{i},num2str(alldata3{j,1}));
            if isequal(alldata3{j,2}, 'default')
                Channels = 1:28;
            else
                Channels = str2num(alldata3{j,2});
            end
            Threshold = alldata3{j,3};
            TimeWindowMinimum = alldata3{j,4};
            TimeWindowMaximum = alldata3{j,5};
            WindowSize = alldata3{j,6};
            WindowStep = alldata3{j,7};
        end
    end

    EEG  = pop_artmwppth( EEG , 'Channel',  Channels, 'Flag', [1 3], 'Threshold', Threshold, 'Twindow', [TimeWindowMinimum  TimeWindowMaximum], 'Windowsize', WindowSize, 'Windowstep', WindowStep ); 
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch_interp_SVT_MW1'], 'gui', 'off'); 

   if ~strcmp(COMP{c}, 'MMN') && ~strcmp(COMP{c}, 'ERN') && ~strcmp(COMP{c}, 'LRP') % not for MMN component
    %Identify segments of EEG with blink artifacts during the stimulus presentation window using the moving window peak-to-peak algorithm with the parameters in the Excel file for this subject
    DimensionsOfFile4 = size(alldata4);
    for j = 1:DimensionsOfFile4(1)
        if isequal(SUB{i},num2str(alldata4{j,1}));
            Channel = alldata4{j,2};
            Threshold = alldata4{j,3};
            TimeWindowMinimum = alldata4{j,4};
            TimeWindowMaximum = alldata4{j,5};
            WindowSize = alldata4{j,6};
            WindowStep = alldata4{j,7};
        end
    end

    EEG  = pop_artmwppth( EEG , 'Channel',  Channel, 'Flag', [1 4], 'Threshold', Threshold, 'Twindow', [TimeWindowMinimum  TimeWindowMaximum], 'Windowsize', WindowSize, 'Windowstep', WindowStep ); 
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch_interp_SVT_MW1_MW2'], 'gui', 'off'); 
   end 

   if ~strcmp(COMP{c}, 'MMN') && ~strcmp(COMP{c}, 'ERN') && ~strcmp(COMP{c}, 'LRP') % not for these components
    %Identify segments of EEG with horizontal eye movement artifacts during the stimulus presentation window using the step like algorithm with the parameters in the Excel file for this subject
    DimensionsOfFile5 = size(alldata5);
    for j = 1:DimensionsOfFile5(1)
        if isequal(SUB{i},num2str(alldata5{j,1}));
            Channel = alldata5{j,2};
            Threshold = alldata5{j,3};
            TimeWindowMinimum = alldata5{j,4};
            TimeWindowMaximum = alldata5{j,5};
            WindowSize = alldata5{j,6};
            WindowStep = alldata5{j,7};
        end
    end

    EEG  = pop_artstep( EEG , 'Channel', Channel, 'Flag', [1 5], 'Threshold',  Threshold, 'Twindow', [TimeWindowMinimum  TimeWindowMaximum], 'Windowsize',  WindowSize, 'Windowstep', WindowStep );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch_interp_SVT_MW1_MW2_SL1'], 'gui', 'off'); 
   end
   
    %Identify segments of EEG with any uncorrected horizontal eye movement artifacts using the step like algorithm with the parameters in the Excel file for this subject
    DimensionsOfFile6 = size(alldata6);
    for j = 1:DimensionsOfFile6(1)
        if isequal(SUB{i},num2str(alldata6{j,1}));
            Channel = alldata6{j,2};
            Threshold = alldata6{j,3};
            TimeWindowMinimum = alldata6{j,4};
            TimeWindowMaximum = alldata6{j,5};
            WindowSize = alldata6{j,6};
            WindowStep = alldata6{j,7};
        end
    end

    EEG  = pop_artstep( EEG , 'Channel', Channel, 'Flag', [1 6], 'Threshold',  Threshold, 'Twindow', [TimeWindowMinimum  TimeWindowMaximum], 'Windowsize',  WindowSize, 'Windowstep', WindowStep );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7, 'setname', [SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch_interp_SVT_MW1_MW2_SL1_SL2'], 'savenew', [Subject_Path SUB{i} '_' COMP{c} '_shifted_ds_reref_ucbip_hpfilt_ica_corr_cbip_elist_bins_epoch_interp_ar.set'], 'gui', 'off'); 

% Additionally save the final preprocessed EEG files to a seperate folder
Preprocessed_Data_Path = fullfile(preprocessed_dataDIR, COMP{c});
mkdir(Preprocessed_Data_Path)
pop_saveset(EEG, ...
    'filename', [SUB{i} '_prepro.set'], ...
    'filepath', Preprocessed_Data_Path);

fprintf('Now finished preprocessing subject  %s\n', SUB{i});
%End subject loop
end

fprintf('Now finished preprocessing all subjects of ERP component: %s\n', COMP{c});

%End ERP component loop
end

%**********************************************************************************************************************************************************************
