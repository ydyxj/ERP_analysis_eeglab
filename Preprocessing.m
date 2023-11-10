%% Start
% Clear memory and the command window.
clear
clc
%  You might need to run EEGLAB and to get some paths in place before you run the script

%% Script Settings
% Set 1 to enable, 0 to disable

clean_up      = 0;    % Cleans up all previous analysis files

import_data   = 0;    % Import raw data files and resample

run_ica       = 0;    % Run ICA, manual input for component rejection is needed after ICA

remove_ica    = 0;    % After run ICA, inspect data,d determine blink components to remove enter it in the list ICA_to_remove

process_data  = 1;    % Filtering, event list, binlister, epoching, artifact rejection, average, lp filter, diff waves

grand_avg     = 1;    % Calculate grant average and grand plots

save_files    = 1;    % Save files for each step. Essential files are saved regardless

plot_figures  = 0;    % Set to 1 to plot figures for each subject


%% Channels to Interpolate
% Plot the data and determine which channels are noisy and needs to be interpolated for each subject
%chan_interpolate = {{},{},{},{'ch10'},{}};
%chan_interpolate = {{'ch10'}};
chan_interpolate = {'Fp1' 'Fp2' 'F7' 'F3' 'Fz' 'F4' 'F8' 'FC5' 'FC1' 'FC2' 'FC6' 'T7' 'C3' 'Cz' 'C4' 'T8' 'CP5' 'CP1' 'CP2' 'CP6' 'P7' 'P3' 'Pz' 'P4' 'P8' 'O1' 'Oz' 'O2'};

%% Global Variables
% Path to the parent folder, which contains the data folders for all subjects
%DIR = /Users/KappenmanLab/ERP_CORE/N170
DIR = fileparts(fileparts(mfilename('fullpath'))); 

%Current_File_Path = /Users/KappenmanLab/ERP_CORE/N170/EEG_ERP_Processing
Current_File_Path = fileparts(mfilename('fullpath'));

ALLERP = buildERPstruct([]); % Initialize the ALLERP structure and CURRENTERP
CURRENTERP = 0;
nsubj = length(subjects); % number of subjects
figScale = [ -300.0 600.0   -300:20:600 ]; % Figure intervals

%% Import Data


subjects = {'S1' 'S2' 'S3' ...
            'S4' 'S5'}; % the variable subjects contains the list of 
                          % subject-specific directories

ALLERP = buildERPstruct([]); % Initialize the ALLERP structure and CURRENTERP
CURRENTERP = 0;
nsubj = numel(subjects); % number of subjects

for s=1:nsubj % Loop through all subjects
        

        % Define subject path based on study directory and subject ID of current subject
        data_path  =  [DIR filesep subjects{s} filesep];

        % import raw data
        %Load the raw continuous EEG data file in .set EEGLAB file format
        EEG = pop_loadset( 'filename', [subjects{s} '_EN.set'], 'filepath', data_path);
        
        EEG  = pop_basicfilter( EEG,  1:33 , 'Boundary', 'boundary', 'Cutoff',  0.1, 'Design', 'butter', 'Filter', 'highpass', 'Order',  2, 'RemoveDC', 'on' );
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1, 'setname', [subjects{s} '_EN'], 'gui', 'off'); 

end

%Rereference to the average of all 33 EEG channels; create a bipolar HEOG channel (HEOG_left minus HEOG_right) and a bipolar VEOG channel (VEOG_lower minus FP2)
EEG = pop_eegchanoperator( EEG, [Current_File_Path filesep 'Rereference_Add_Uncorrected_Bipolars_N170.txt']);




        %re-reference to average
        EEG=pop_chanedit(EEG, 'lookup','C:\\eeglab14_1_2b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');
        %[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

        EEG = pop_reref( EEG, [],'refloc',struct('labels',{'Cz'},'type',{''},'theta',{0},'radius',{0},'X',{5.2047e-15},'Y',{0},'Z',{85},'sph_theta',{0},'sph_phi',{90},'sph_radius',{85},'urchan',{32},'ref',{''},'datachan',{0}));
        EEG.setname=[EEG.setname  '_reref'];
        if (save_files)
            EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', data_path);
        end

        % HP filter at 0.1
        EEG  = pop_basicfilter( EEG,  1:32 , 'Boundary', 'boundary', 'Cutoff', 0.1, 'Design', 'butter', 'Filter', 'highpass', 'Order',  2, 'RemoveDC', 'on' );
        EEG.setname = [EEG.setname '_hpfilt'];
        EEG = pop_saveset(EEG, 'filename', [EEG.setname '.set'], 'filepath', data_path);

    end
end


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%load EEG data 
EEG = pop_loadbv('C:\Personal\Yan\eeg data\Priming\sub2\En\', 'Priming_sub0002.vhdr', [1 1948580], [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31]);
%set name for dataset
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','Sub2 EN','gui','off'); 

%edit channel location
EEG=pop_chanedit(EEG, 'lookup','C:\\eeglab14_1_2b\\plugins\\dipfit2.3\\standard_BESA\\standard-10-5-cap385.elp');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%view channel data(scroll)
pop_eegplot( EEG, 1, 0, 1);

%Create EVENTLIST
EEG  = pop_creabasiceventlist( EEG , 'AlphanumericCleaning', 'on', 'BoundaryNumeric', { -99 }, 'BoundaryString', { 'boundary' }, 'Eventlist', 'C:\Personal\Yan\eeg data\Priming\sub2\En\elist.txt' ); 
 

%Assign  Binlister
EEG  = pop_binlister( EEG , 'BDF', 'C:\Personal\Yan\eeg data\Priming\binlister_demo_1.txt', 'ExportEL', 'C:\Personal\Yan\eeg data\Priming\sub5\En\elist 1.txt', 'IndexEL',  1, 'SendEL2', 'EEG&Text', 'Voutput', 'EEG' ); 
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%Epoch
EEG = pop_epochbin( EEG , [-300.0  600.0],  'pre'); 
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'setname','Sub2 EN_test','savenew','C:\\Personal\\Yan\\eeg data\\Priming\\sub2\\En\\Sub2 EN_test','gui','off'); 


EEG  = pop_artmwppth( EEG , 'Channel',  4:31, 'Flag',  1, 'Threshold',  100, 'Twindow', [ -200 299], 'Windowsize',  200, 'Windowstep',  100 ); % GUI: 01-Dec-2021 14:27:19
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname','Sub2 EN_elist_bins_epoch_artifical','savenew','C:\\Personal\\Yan\\eeg data\\Priming\\sub2\\En\\Sub2 EN_elist_bins_epoch_artifical.set','gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_runica(EEG, 'extended',1,'interupt','on');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
pop_selectcomps(EEG, [1:31] );
EEG  = pop_artmwppth( EEG , 'Channel',  6:31, 'Flag',  1, 'Threshold',  100, 'Twindow', [ -200 299], 'Windowsize',  200, 'Windowstep',  100 ); % GUI: 01-Dec-2021 14:35:26
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'gui','off'); 
EEG  = pop_artmwppth( EEG , 'Channel',  6:31, 'Flag',  1, 'Threshold',  100, 'Twindow', [ -200 299], 'Windowsize',  200, 'Windowstep',  100 ); % GUI: 01-Dec-2021 14:37:28
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'gui','off'); 
EEG  = pop_basicfilter( EEG,  1:31 , 'Cutoff',  30, 'Design', 'butter', 'Filter', 'lowpass', 'Order',  6 ); % GUI: 01-Dec-2021 14:39:39
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7,'savenew','C:\\Personal\\Yan\\eeg data\\Priming\\sub2\\En\\Sub2 EN_elist_bins_epoch_artifical_ar_filt.set','gui','off'); 
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
EEG = pop_rejepoch( EEG, [3:6 8 9 11 12:16 18 19 21 27 28 30 31 32 34 37 39 40 43 45 46 47 50 51 53 57 58 59 61 62:70 72 73:75 77 78:80 84 87 88 91 92 96 97:100 107 110 112 116 117 119 120 121 127 132 133 136 140 142:2:146 156 158 162 163 165 167 171 172 176 177:180 182 183:193 195 196 198 199 201 202 203 206 207 211 212:214 216 217:219 223 225 227 228 229 231 232 234 235:250 252 256 257 258 260 267 277 278 279 286 291 292 294] ,0);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 8,'setname','Sub2 EN_elist_bins_epoch_artifical_reject','savenew','C:\\Personal\\Yan\\eeg data\\Priming\\sub2\\En\\Sub2 EN_elist_bins_epoch_artifical_reject.set','gui','off'); 
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 9,'retrieve',9,'study',0); 
STUDY = []; CURRENTSTUDY = 0; ALLEEG = []; EEG=[]; CURRENTSET=[];