%% 1. Setup
clear
close all


cd 'E:\LaCie\DATA_KD\161005_E'
 
%load WorkspaceForProcessing_180705_evp001_test.mat
load Workspace_161005_E_evp001_2018-07-25.mat
    % Available variables, 'EV.tp' 'Filename' 'LFP' 'MUA' and saveName'
    
%% 2. Trigger windows

[stimWin] = triggerStim(EV.tp,LFP);

%% 3. Split the electrodes
el.C.full   = stimWin.stimLFP_full(:,:,1:24);
el.C.bl     = stimWin.stimLFP_bl(:,:,1:24);
el.C.target = stimWin.stimLFP_target(:,:,1:24);
el.D.full   = stimWin.stimLFP_full(:,:,25:48);
el.D.bl     = stimWin.stimLFP_bl(:,:,25:48);
el.D.target = stimWin.stimLFP_target(:,:,25:48);

%Baseline correct raw LFP
%elC
avgBl.C.txchan = squeeze(median(el.C.bl,1)); % average across trials. produce time by channel
avgBl.C.chan = median(avgBl.C.txchan,1); % average across time
avgBl.C.final = median(avgBl.C.chan,2); % average across channels
el.C.blCorrect = el.C.full-avgBl.C.final; %subtract from every value in full

% elD
avgBl.D.txchan = squeeze(median(el.D.bl,1));
avgBl.D.chan = median(avgBl.D.txchan,1);
avgBl.D.final = median(avgBl.D.chan,2);
el.D.blCorrect = el.D.full-avgBl.D.final;

%average across trials after baselin correcting
 viewLFP.C = squeeze(median(el.C.blCorrect,1));
 viewLFP.D = squeeze(median(el.D.blCorrect,1));

%time on the x axis el contact on y axis
figure
spLFP_elC = subplot(1,2,1);
    plot(viewLFP.C)
    title('Baseline correct LFP elC');
spLFP_elD = subplot(1,2,2);
    plot(viewLFP.D)
    title('Baseline correctLFP elD');
    
   
    
%%
%CONGRATS... it finished
load gong
sound(y,Fs)
