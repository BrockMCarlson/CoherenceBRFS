%181113 current problems getting the EV.A-D tp and the LFP data to work in
%this coherence script. Problems with syntax in structure variables. 

% Current error at line 37: "Undefined function 'sort' for input arguments 
% of type 'struct'."



%% 1. Setup
clear
close all

%Variable inputs
sink = 16;

cd 'E:\LaCie\DATA_KD\161005_E'
 
%load WorkspaceForProcessing_180705_evp001_test.mat
load BRFSWorkspace_161005_E_brfs001ns2_2018-11-13.mat
    % Available variables, 'EV.tp' 'Filename' 'LFP' 'MUA' and saveName'
 
%% 2. Trigger windows

[stimWin.A] = triggerStim(EV.A,LFP);
[stimWin.B] = triggerStim(EV.B,LFP);
[stimWin.C] = triggerStim(EV.C,LFP);
[stimWin.D] = triggerStim(EV.D,LFP);

%% 3. Split the electrodes

el.C.full   = stimWin.A.stimLFP_full(:,:,1:24);
el.C.bl     = stimWin.A.stimLFP_bl(:,:,1:24);
el.C.target = stimWin.A.stimLFP_target(:,:,1:24);
% % % % % el.D.full   = stimWin.stimLFP_full(:,:,25:48);
% % % % % el.D.bl     = stimWin.stimLFP_bl(:,:,25:48);
% % % % % el.D.target = stimWin.stimLFP_target(:,:,25:48);
 
%% 4. Calculate the Coherence
[coher,f] = multiCoher(sink,el);

cxy_target  = coher.cxy.target;
cxy_bl      = coher.cxy.bl;

% 5. Average coherence across trials
cxy_target_avg = median(cxy_target,1);
cxy_bl_avg     = median(cxy_bl,1);

% 6.  Baseline subtracted coherence
% subtract coherence baseline averaged across trials from the coherence of
%each trial of the target window data

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

  


      
%Baseline correct the LFP????
 
 
% 10. Plot 
close all

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

