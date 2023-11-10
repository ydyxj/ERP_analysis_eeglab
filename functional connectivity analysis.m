 

%% Preprocessing

clear all; close all; clc

%% 指定相关信息
Data_Dir = uigetdir([],'Path of the EEG datasets'); %指定待分析的脑电数据所在的路径；可以是预处理后的分段脑电数据。
                                                    % 也可以是CSD转换后的分段数据；
                                                    % 建议一个条件或一个组别一个文件夹，每个文件夹分别运行该代码

Output_Dir = uigetdir([],'Path to store the measures'); %计算得到的指标文件保存的路径

band = inputdlg('the limits of band'); %指定所要分析的频率的范围（单位是Hz）
band = str2num(band{1}); %将band变量由字符转换为数值；注意：运行本脚本时MATLAB的路径中不应该有HERMES工具箱

bandname = inputdlg('the name of the band you computed'); % 指定前述待分析的频段的名称
bandname = bandname{1};

prefix = inputdlg('the prefix of computed measures');%指定保存的指标文件的前缀，建议是条件名或组别名
prefix = prefix{1};

%% 获取EEG数据路径中包含的set文件
Dir_Data = dir(fullfile(Data_Dir,'*.set')); 
FileNames = {Dir_Data.name};

%% compute 1st measure (coherence)
for subj = 1:numel(FileNames)  
    EEG = pop_loadset('filename',FileNames{1,subj},'filepath',Data_Dir); % 载入某个被试的数据
    EEG = eeg_checkset( EEG );
    N = EEG.pnts; %每段的长度（点数）
    SampleRate = EEG.srate; %取样率 sampling rate
    NFFT = 2^nextpow2(N); %大于每段长度的、最小的2的N次方
    Freq = SampleRate/2*linspace(0,1,NFFT/2+1); %频率轴 
    for chan = 1:size(EEG.data,1)
        for epochs = 1:size(EEG.data,3)
            ffts(:,chan,epochs) = fft(hanning(N).*squeeze(EEG.data(chan,:,epochs))', NFFT);% 对该被试每个电极的每个分段进行FFT            
        end
    end
    for x = 1:size(EEG.data,1)
            for y = 1:size(EEG.data,1)
                fx = squeeze(ffts(:,x,:));
                Pxx = fx.*conj(fx)/N;
                MeanPx = mean(Pxx,2); % 计算coherence时，x电极的power
                fy = squeeze(ffts(:,y,:));
                Pyy = fy.*conj(fy)/N; 
                MeanPy = mean(Pyy,2); % 计算coherence时，y电极的power
                Pxy = fx.*conj(fy)/N;
                MeanPxy = mean(Pxy,2); %% Sxy，上述两个电极的交叉谱
                C = (abs(MeanPxy).^2)./(MeanPx.*MeanPy); % 相干
                coh(:,x,y,subj) = C; % coherence，频率*电极*电极*被试
            end
    end

    clear ffts
end

%-------------------------------------------------------------------------
coh = coh(1:NFFT/2 + 1,:,:,:);  

idx = dsearchn(Freq', band'); %确定频段上下限在Freq即频率轴中的位置
coh = squeeze(mean(coh(idx(1,1):idx(2,1),:,:,:),1));%计算某个频段的平均coh

save(strcat(Output_Dir,'\',prefix,'_',bandname,'_coh.mat'),'coh');

%% 2nd and 3rd measures (phase-locking value and phase lag index)   
for subj = 1:numel(FileNames)
    EEG = pop_loadset('filename',FileNames{1,subj},'filepath',Data_Dir);
    EEG = eeg_checkset( EEG );
    eeg_filtered = eegfilt(reshape(EEG.data, [size(EEG.data,1) size(EEG.data,2)*size(EEG.data,3)]),...
                   EEG.srate,band(1,1),band(1,2),0,3*fix(EEG.srate/band(1,1)),0,'fir1',0); % 对载入的数据进行带通滤波；使用eegfilt进行带通滤波时需要将分段数据重新变为连续数据 
         
    for channels = 1:size(EEG.data,1)
        band_phase(channels,:) = angle(hilbert(eeg_filtered(channels,:))); %逐个分段进行Hilbert变换，并提取相位
    end    
    perc10w =  floor(size(band_phase,2)*0.1);% 确定数据长度10%是多少个样本点
    band_phase = band_phase(:,perc10w+1:end-perc10w); %因Hilbert变换对数据首尾相位估算不准确，顾去掉前10%和后10%样本点的相位
    epoch_num = floor(size(band_phase,2)/size(EEG.data,2)); % 确定剩余的样本点如果转换为分段数据，可以分成多少段
    band_phase = band_phase(:,1:epoch_num*size(EEG.data,2)); % 依据可以分成的段数，截取数据
    band_phase = reshape(band_phase,[size(EEG.data,1) size(EEG.data,2) epoch_num]);% 将数据重新转换为三维：电极*样本点*分段
    
    for x = 1:size(band_phase,1)
         for y = 1:size(band_phase,1)
             for epochs = 1:size(band_phase,3)
                 x_phase = squeeze(band_phase(x,:,epochs)); % 提取电极对中第一个电极在某个分段的相位
                 y_phase = squeeze(band_phase(y,:,epochs)); % 提取电极对中第二个电极在某个分段的相位
                 rp = x_phase - y_phase; % 计算两个电极在某个分段的相位差
                 %%% PLV
                 sub_plv(x,y,epochs) = abs(sum(exp(1i*rp))/length(rp)); % 计算某个被试某个电极对在某个分段的PLV
                 %%% PLI
                 sub_pli(x,y,epochs) = abs(mean(sign((abs(rp)- pi).*rp))); % 计算某个被试某个电极对在某个分段的PLI
             end
         end
    end
   pli(:,:,subj) = mean(sub_pli,3); % 对该被试各个分段的PLI计算平均值；pli变量维度是电极*电极*被试
   plv(:,:,subj) = mean(sub_plv,3); % 对该被试各个分段的PLV计算平均值；plv变量维度是电极*电极*被试
   clear band_phase sub_pli sub_plv
end

save(strcat(Output_Dir,'\',prefix,'_',bandname,'_plv.mat'),'plv');
save(strcat(Output_Dir,'\',prefix,'_',bandname,'_pli.mat'),'pli');

%% compute 4th measure (weighted phase lag index)
for subj = 1:numel(FileNames)
    EEG = pop_loadset('filename',FileNames{1,subj},'filepath',Data_Dir);
    EEG = eeg_checkset( EEG );
    eeg_filtered = eegfilt(reshape(EEG.data, [size(EEG.data,1) size(EEG.data,2)*size(EEG.data,3)]),...
                   EEG.srate,band(1,1),band(1,2),0,3*fix(EEG.srate/band(1,1)),0,'fir1',0); 
    
    for channels = 1:size(EEG.data,1)
        band_hilbert(channels,:) = hilbert(eeg_filtered(channels,:)); % 对该被试每个通道的信号进行hilbert变换，得到解析信号（a + b*i）；计算wPLI不需要提取相位
    end 
    perc10w =  floor(size(band_hilbert,2)*0.1);
    band_hilbert = band_hilbert(:,perc10w+1:end-perc10w);
    epoch_num = floor(size(band_hilbert,2)/size(EEG.data,2));
    band_hilbert = band_hilbert(:,1:epoch_num*size(EEG.data,2));
    band_hilbert = reshape(band_hilbert,[size(EEG.data,1) size(EEG.data,2) epoch_num]); 
    
    for x = 1:size(band_hilbert,1)
         for y = 1:size(band_hilbert,1)
             for epochs = 1:size(band_hilbert,3)
                 x_hilbert = band_hilbert(x,:,epochs);
                 y_hilbert = band_hilbert(y,:,epochs);
                 crossspec = x_hilbert.* conj(y_hilbert); % 交叉谱
                 crossspec_imag = imag(crossspec); % 交叉谱的虚部
                 sub_wpli(x,y,epochs) = abs(mean(crossspec_imag))/mean(abs(crossspec_imag)); % 计算某个被试某个电极对在某个分段的wPLI
             end
         end
    end
    wpli(:,:,subj) = mean(sub_wpli,3); 
    clear band_hilbert sub_wpli
end
wpli(isnan(wpli)) = 0;% 因为对角线上的wpli数值为nan，故进行这个操作
save(strcat(Output_Dir,'\',prefix,'_',bandname,'_wpli.mat'),'wpli');  
%% plot  所有被试 
figure;
for i=1:size(plv,3)
    subj_fc = plv(:,:,i);
    subplot(3,5,i);
    imagesc(subj_fc);
    title(['subj',num2str(i)]);
end

%% 组平均
group1_fc = mean(plv,3);
figure;
imagesc(group1_fc);
title('Group-Average  FC');
