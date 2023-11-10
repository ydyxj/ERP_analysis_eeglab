clear all; clc; close all

%% part1: compute group-level ERP

Subj = [1:10]; %% subject numbers
for i = 1:length(Subj)
    setname = strcat('s',num2str(i),'EN','.set'); %% filename of set file
    setpath = 'D:\MyWorkSpace\Matlab\SiYingPeiXun\25EEG_day2\Example_data\';

    EEG = pop_loadset('filename',setname,'filepath',setpath); %% load the data
    EEG = eeg_checkset( EEG );


    % �?总每个被试个体水平�?�?平�?好�?数据
    % 10 59 3000  sub chl time
    EEG_avg(i,:,:) = squeeze(mean(EEG.data,3)); %% single-subject ERPs; EEG_avg dimension: subj*channel*time
end

save('Group_level_ERP.mat','EEG_avg');  %% save the data of subjects

%--------------------------------------------------------------------------

%% part2: plot group-level ERP

Cz = 17; % select the channel to plot (display maximum response)
figure;
% 提取�?有被试Cz电极所有时间点�?数据
% mean_data 是Cz电极所有时间点�?�?水平平�?数据
mean_data = squeeze(mean(EEG_avg(:,Cz,:),1)); %% select the data at Cz, average across subjects, mean_data: 1*3000

plot(EEG.times, mean_data,'k','linewidth', 1.5); %% plot the waveforms
set(gca,'YDir','reverse'); %% reverse the direction of Y axis

axis([-350 700 -15 10]);  %% define the region to display
% xlim([-350 700]);  %% define the region of X axis
% ylim([-15 10]); %% define the region of Y axis

title('Group-level at Cz','fontsize',16); %% specify the figure name
xlabel('Latency (ms)','fontsize',16); %% name of X axis
ylabel('Amplitude (uV)','fontsize',16);  %% name of Y axis


%--------------------------------------------------------------------------
%% part3: plot the scalp maps at dominant peaks

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

figure; %subplot(2,2,1) one row and two columns worth of figures. p=1 means that you wish to place the plot in the left most column
subplot(221); topoplot(mean(P1_amplitude),EEG.chanlocs,'maplimits',[-15 15]); title('N2 Amplitude','fontsize',16); %% P1 scalp map (group-level)
subplot(222); topoplot(mean(N1_amplitude),EEG.chanlocs,'maplimits',[-15 15]); title('P2 Amplitude','fontsize',16); %% N1 scalp map (group-level
subplot(223); topoplot(mean(N2_amplitude),EEG.chanlocs,'maplimits',[-15 15]); title('N2 Amplitude','fontsize',16); %% N2 scalp map (group-level)
subplot(224); topoplot(mean(P2_amplitude),EEG.chanlocs,'maplimits',[-15 15]); title('P2 Amplitude','fontsize',16); %% P2 scalp map (group-level)


%%------------------------------------------------------------------------
%% part4: series of scalp mps

time_interval = [0:50:450]; %% specify the time intervals to display (to be changed)
figure; 
for i = 1:length(time_interval)
    latency_range = [time_interval(i) time_interval(i)+50]; %% lower and upper limits
    latency_idx = find((EEG.times>=latency_range(1))&(EEG.times<=latency_range(2))); %% interval of the specific regions
    Amplitude = squeeze(mean(mean(EEG_avg(:,:,latency_idx),1),3)); %% 1*channel (averaged across subjects and interval)
    subplot(3,3,i); 
    topoplot(Amplitude,EEG.chanlocs,'maplimits',[-10 10]); %% topoplot(Amplitude,EEG.chanlocs);
    setname = strcat(num2str(latency_range(1)),'--',num2str(latency_range(2)),'ms'); %% specify the name of subplots
    title(setname,'fontsize',16); %% display the names of subplots
end

