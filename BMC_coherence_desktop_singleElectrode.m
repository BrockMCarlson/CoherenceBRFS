%% 1. Setup
clear
close all

%Variable inputs
sink = 16;
top = sink-2;
bottom = sink+2;

cd 'F:\LaCie\DATA_KD\161005_E'
 
%load WorkspaceForProcessing_180705_evp001_test.mat
load Workspace_161005_E_evp001_2018-07-25.mat
% Available variables, 'EV.tp' 'Filename' 'LFP' 'MUA' and saveName'
 
 
%% 2. Trigger windows
clear stimLFP
pre = 256;   % pre-stim time (baseline) in ms
pre_bl = 256;
post = 612; % post-stim time in ms
postTargetStart = 100;
postTargetStop = postTargetStart+512;
 
% TRIGGER LFP TO STIM ON
for tr = 1:length(EV.tp) % trigger to stim-on times for all trials
    stimtm = round(EV.tp(tr,1)/30) ;% divide by 30 to convert to 1kHz. Note, LFP and MUA already in 1kHZ
    refwin = stimtm-pre:stimtm+post;
    refwin_bl =  stimtm-pre_bl:stimtm;
    refwin_target = stimtm+postTargetStart:stimtm+postTargetStop;
    stimLFP(tr,:,:)        = LFP(refwin,:);
    stimLFP_bl(tr,:,:)     = LFP(refwin_bl,:);
    stimLFP_target(tr,:,:) = LFP(refwin_target,:);
    
end


%% 3. Split the electrodes
LFP_elC        = stimLFP(:,:,1:24);
LFP_elC_bl     = stimLFP_bl(:,:,1:24);
LFP_elC_target = stimLFP_target(:,:,1:24);
LFP_elD        = stimLFP(:,:,25:48);
LFP_elD_bl     = stimLFP_bl(:,:,25:48);
LFP_elD_target = stimLFP_target(:,:,25:48);
 
 
%% 4. Calculate the Coherence
clear x_bl x_target y_bl b_target 

%parameters
window=128;
noverlap=127;
nfft=128;
fs=1000;
f = 1:500;

% elC
% x and y comparision setting
% ********* Change to proper electrode and channel*********
clear x_bl x_target y_bl b_target 
x_target = LFP_elD_target(:,:,sink);
x_bl    = LFP_elD_bl(:,:,sink);
y_target = LFP_elD_target(:,:,top);
y_bl     = LFP_elD_bl(:,:,top);
cxy_target = nan(size(x_target,1),65);
cxy_bl    = nan(size(x_bl,1),65);


%run mscohere
clear i 
for i = 1:size(x_bl,1)
    [cxy_target(i,:),f] = mscohere(x_target(i,:),y_target(i,:),window,noverlap,nfft,fs);
    [cxy_bl(i,:),f]     = mscohere(x_bl(i,:),y_bl(i,:),window,noverlap,nfft,fs);
end

% 5. Average coherence across trials
cxy_target_avg = median(cxy_target,1);
cxy_bl_avg     = median(cxy_bl,1);

% 6.  Baseline subtracted coherence
% subtract coherence baseline averaged across trials from the coherence of
%each trial of the target window data
for i = 1:size(cxy_target,1) %the coherence of every trial
       cxy_target_blCorrect(i,:) = cxy_target(i,:)-cxy_bl_avg; %419x65(i,:) - 1x65
end
% 7. Average coherence across trials of baseline subtracted data
cxy_targetBlCorrect_avg = median(cxy_target_blCorrect,1);

% 8. Find the confidence interval for coherence across trials
% Note, this needs to be done before data is averaged across trials.
% However, the reason that it appears here in the code is that the baseline
% corrected coherence for the target window needs the baseline data first to be
% averaged across trials.
reps = 5000;
 % ***CALLS FROM UNAVERAGED DATA***
cxy_targetBlCorrect_CI = bootci(reps,@median,cxy_target_blCorrect);
cxy_target_CI = bootci(reps,@median,cxy_target); % no baseline correct
cxy_bl_CI     = bootci(reps,@median,cxy_bl);

  


      
% 9. Baseline correct raw LFP and view LFP
%Baseline correct raw LFP

%elC
baseline_txchan_elC = squeeze(median(LFP_elC_bl,1));
baseline_chan_elC = median(baseline_txchan_elC,1);
baseline_elC = median(baseline_chan_elC,2);
stimLFP_blCorrect_elC = LFP_elC-baseline_elC;

% elD
baseline_txchan_elD = squeeze(median(LFP_elD_bl,1));
baseline_chan_elD = median(baseline_txchan_elD,1);
baseline_elD = median(baseline_chan_elD,2);
stimLFP_blCorrect_elD = LFP_elD-baseline_elD;

%Look at raw LFP
 viewLFP_elC = squeeze(median(stimLFP_blCorrect_elC,1));
 viewLFP_elD = squeeze(median(stimLFP_blCorrect_elD,1));


 
 
% 10. Plot 
close all

% A. Plot raw LFP
figure(1)
spLFP_elC = subplot(1,2,1);
    plot(viewLFP_elC)
    title('Baseline correct LFP elC');
spLFP_elD = subplot(1,2,2);
    plot(viewLFP_elD)
    title('Baseline correctLFP elD');

% B. Plot coherence tests
figure(2)
trgt_blAvgSub = subplot(3,1,1);
    plot(f,cxy_targetBlCorrect_avg); hold on;
    plot(f,cxy_targetBlCorrect_CI(1,:),'linestyle','--'); hold on;
    plot(f,cxy_targetBlCorrect_CI(2,:),'linestyle',':'); hold on;
    xlim([0 100]);
    title('elD. Sink to bottom. Target, baseline subed');hold off;

trgt = subplot(3,1,2);
    plot(f,cxy_target_avg);  hold on;
    plot(f,cxy_target_CI(1,:),'linestyle','--'); hold on;
    plot(f,cxy_target_CI(2,:),'linestyle',':'); hold on;
    xlim([0 100]);
    title('elD. Sink to bottom. Target'); hold off;

bl = subplot(3,1,3);
    plot(f,cxy_bl_avg); hold on;
    plot(f,cxy_bl_CI(1,:),'linestyle','--'); hold on;
    plot(f,cxy_bl_CI(2,:),'linestyle',':'); hold on;
    xlim([0 100]);
    title('elD. Sink to bottom. Baseline');hold off;

%%
%CONGRATS... it finished
load gong
sound(y,Fs)

    
