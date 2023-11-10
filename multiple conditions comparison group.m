%%-------------------------------------------------------------------------
%Condition comparison

clear all; clc; close all 

%%  part1: plot the waveforms for different conditions

Cond = {'congruent','baseline','incongruent'}; %% condition name

for i = 1:length(Subj)
     setname = strcat('s',num2str(i),'EN','.set'); %% name of the set file
     setpath = 'D:\Study\eeg data\audio priming\rawdata'; %% filepath of the set file
     EEG = pop_loadset('filename',setname,'filepath',setpath);  %% load the data into EEG
     EEG = eeg_checkset( EEG );
    
    for j = 1:length(Cond)

        EEG_new = pop_epoch( EEG, Cond(j), [-0.6  0.8], 'newname', 'Merged datasets', 'epochinfo', 'yes'); %% epoch by conditions, input to EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_new = pop_rmbase( EEG_new, [-600  0]); %% baseline correction for EEG_new
        EEG_new = eeg_checkset( EEG_new );
        EEG_avg(i,j,:,:) = squeeze(mean(EEG_new.data,3));  %% average across trials for EEG_new, EEG_avg dimension: subj*cond*channel*time
    end 

Cz = 17; %% channel to display
mean_data = squeeze(mean(EEG_avg(:,:,Cz,:),1)); %% select data at Cz, and average across subjects, mean_data: cond*times
figure; 
plot(EEG.times, mean_data,'linewidth', 1.5); %% plot waveforms for different conditions
set(gca,'YDir','reverse');  %% reverse Y axis
axis([-500 1000 -35 25]);  %% define the region to display
title('Group level data','fontsize',16); 
xlabel('Latency (ms)','fontsize',16);
ylabel('Amplitude (uV)','fontsize',16); 
legend(Cond)

%--------------------------------------------------------------------------
%% scalp maps of dominant peak for different conditions

P1_peak = 50;
N1_peak = 100; 
N2_peak = 207; 
P2_peak = 374; %% dominant peaks on waveforms

P1_interval = find((EEG.times>=197)&(EEG.times<=217)); %% define the P1 intervals [peak-10 peak+10]
N1_interval = find((EEG.times>=364)&(EEG.times<=384)); %% define the N1 intervals [peak-10 peak+10]
N2_interval = find((EEG.times>=197)&(EEG.times<=217)); %% define the N2 intervals [peak-10 peak+10]
P2_interval = find((EEG.times>=364)&(EEG.times<=384)); %% define the P2 intervals [peak-10 peak+10]

P1_amplitude = squeeze(mean(EEG_avg(:,:,P1_interval),3));  %% P1 amplitude for each subject and each channel
N1_amplitude = squeeze(mean(EEG_avg(:,:,N1_interval),3));   %% N1 amplitude for each subject and each channel
N2_amplitude = squeeze(mean(EEG_avg(:,:,N2_interval),3));  %% N2 amplitude for each subject and each channel
P2_amplitude = squeeze(mean(EEG_avg(:,:,P2_interval),3));   %% P2 amplitude for each subject and each channel

figure; %% divide the panel into 2 rows and 3 colums
for i = 1:3
    % 提取P1成�?该条件�?�?有被�?  当前条件  �?有�?�道�?数据
    P1_data = squeeze(mean(P1_amplitude(:,i,:),1)); %% average across subjects
    subplot(2,3,i); 
    topoplot(P1_data,EEG.chanlocs,'maplimits',[-15 15]); 
    colorbar;
    titlename = strcat('P1 Amplitude ',Cond(i))
    title(titlename ,'fontsize',16); %% plot N2 scalp map (group-level)

% 提取N1成�?该条件�?�?有被�?  当前条件  �?有�?�道�?数据
    N1_data = mean(N1_amplitude(:,i,:),1); %% average across subjets
    subplot(2,3,i+3); 
    topoplot(P2_data,EEG.chanlocs,'maplimits',[-15 15]); 
    colorbar;
    title(titlename,'fontsize',16); %% plot P2 scamp map (group-level)
end



%%---------------------------------------------------------------------------
%% 比�?单个个体�??条件�? 时间上的差�?
clear all; clc; close all

Subj = [1:10]; 
Cond = {'congruent','baseline','incongruent'};

%% compute averaged data
for i = 1:length(Subj)
    setname = strcat('s',num2str(i),'EN','.set'); 
    setpath = 'D:\MyWorkSpace\Matlab\SiYingPeiXun\25EEG_day2\Example_data\';
    EEG = pop_loadset('filename',setname,'filepath',setpath); 
    EEG = eeg_checkset( EEG );

    for j = 1:length(Cond)
        EEG_new = pop_epoch( EEG, Cond(j), [-0.6 0.8], 'newname', 'Merged datasets pruned with ICA', 'epochinfo', 'yes'); 
        EEG_new = eeg_checkset( EEG_new );
        EEG_new = pop_rmbase( EEG_new, [-600  0]); 
        EEG_new = eeg_checkset( EEG_new );
        EEG_avg(i,j,:,:)=squeeze(mean(EEG_new.data,3));  %% subj*cond*channel*timepoints
    end 
end


%% T test
%%  point-by-point paried t-test  across multiple time points 

data_test = squeeze(EEG_avg(:,:,17,:)); %% select the data at Cz, data_test: subj*cond*time
for i = 1:size(data_test,3)
    data_1 = squeeze(data_test(:,1,i)); %% select condition congruent for each time point
    data_2 = squeeze(data_test(:,3,i)); %% select condition incongruent for each time point

    % �?验两个条件�?差�?  用配对样本T�?�?  ttest
    % �?验两�?样本之间�?差异，要做独立�?�本T�?�?  ttest2

    [h p] = ttest(data_1,data_2); %% ttest comparison
    P_ttest(i) = p; %% save the p value from ttest


end
figure; 
subplot(211); plot(EEG.times,squeeze(mean(data_test(:,1,:),1)),'b'); %% plot the average waveform for congruent
hold on; plot(EEG.times,squeeze(mean(data_test(:,3,:),1)),'r'); %% plot the average waveform for incongruent
subplot(212); plot(EEG.times,P_ttest); ylim([0 0.05]); %%plot the p values from ttest


%% -------------------------------------------------------------------------

%% compare channel difference

%% point-by-point paried t-test  across multiple channels
% define ROI --> local component
test_idx = find((EEG.times>=197)&(EEG.times<=217)); %% define the intervals
data_test = squeeze(mean(EEG_avg(:,:,:,test_idx),4)); %% select the data in [197 217]ms, subj*cond*channel
for i = 1:size(data_test,3)
    data_1 = squeeze(data_test(:,1,i)); %% select condition L3 for each channel
    data_2 = squeeze(data_test(:,3,i)); %% select condition L4 for each channel
    [h,p,ci,stats] = ttest(data_1,data_2); %% ttest comparison
    P_ttest2(i) = p; %% save the p value from ttest
    T_ttest2(i) = stats.tstat; 
end

figure; 
subplot(141); 
topoplot(squeeze(mean(data_test(:,1,:),1)),EEG.chanlocs,'maplimits',[-20 20]); 
subplot(142); 
topoplot(squeeze(mean(data_test(:,3,:),1)),EEG.chanlocs,'maplimits',[-20 20]); 
subplot(143); 
topoplot(T_ttest2,EEG.chanlocs); 
subplot(144); 
topoplot(P_ttest2,EEG.chanlocs,'maplimits',[0 0.05]); 

%% ------------------------------------------------------------------------

%% ANOVA for 2 more conditions
%% point-by-point repeated measures of ANOVA across time points
data_test = squeeze(EEG_avg(:,:,17,:)); %% select the data at Cz, data_test: subj*cond*time
for i = 1:size(data_test,3)
    data_anova = squeeze(data_test(:,:,i)); %% select the data at time point i
    [p, table] = anova_rm(data_anova,'off');  %% perform repeated measures ANOVA
    P_anova(i) = p(1); %% save the data from ANOVA
end

mean_data = squeeze(mean(data_test,1)); %% dimension: cond*time
figure; 
subplot(211);plot(EEG.times, mean_data,'linewidth', 1.5); %% waveform for different condition 
set(gca,'YDir','reverse');
axis([-500 800 -35 25]);
subplot(212);plot(EEG.times,P_anova); axis([-500 800 0 0.05]); %% plot the p values from ANOVA


%% -----------------------------------------------------------------------

%% point-by-point repeated measures of ANOVA across channels

test_idx = find((EEG.times>=197)&(EEG.times<=217)); %% define the intervals
data_test = squeeze(mean(EEG_avg(:,:,:,test_idx),4)); %% select the data in [197 217]ms, subj*cond*channel
for i = 1:size(data_test,3)
    data_anova = squeeze(data_test(:,:,i)); %% select the data at channel i
    [p, table] = anova_rm(data_anova,'off');  %% perform repeated measures ANOVA
    P_anova2(i) = p(1); %% save the data from ANOVA
    F_anova2(i) = table{2,5};
end
figure; 
for i = 1:4
    subplot(1,5,i); 
    topoplot(squeeze(mean(data_test(:,i,:),1)),EEG.chanlocs,'maplimits',[-20 20]); 
end
% subplot(1,5,5); topoplot( P_anova2,EEG.chanlocs,'maplimits',[0 0.05]); 
subplot(1,5,5); topoplot( F_anova2,EEG.chanlocs); 





