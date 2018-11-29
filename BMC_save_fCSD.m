%BMC plot CSD
%loads a matlab workspace created in ReadInNS6 to use variables for
%processing and then plotting. Available variables, 'EV.tp' 'Filename'
%'LFP' 'MUA' and saveName'
%trigger to stimulus on based on event codes
%average LFP across trials
    %run calcCSD
    % convert uV to nA/mm^3
    % subtract baseline
    % 1/2 wave rectify
    % pad out to same number of channels as LFP
    % filter CSD
    %plot

clear
close all
cd 'E:\LaCie\DATA_KD\161007_E\'

%load WorkspaceForProcessing_180705_evp001_test.mat
load BRFSWorkspace_161007_E_brfs001ns2_2018-11-15.mat
% Available variables, 'Cond' 'EV' 'extension' 'Filename' 'LFP' 'saveName'
% and 'Unq_cond' 


%% 4. PROCESS
clear stimLFP
pre = 100;   % pre-stim time (baseline) in ms
post = 500; % post-stim time in ms

% TRIGGER MUA & LFP TO STIM ON
for tr = 1:length(EV.D) % trigger to stim-on times for all trials
    stimtm = round(EV.D(tr,2)/30) ;% divide by 30 to convert to 1kHz. Note, LFP and MUA already in 1kHZ
    refwin = stimtm-pre:stimtm+post;
    stimLFP(tr,:,:) = LFP(refwin,:);
end

% COMPUTE AVERAGE ACROSS TRIALS
avgLFP = squeeze(mean(stimLFP,1));
avgLFP_elC = avgLFP(:,1:24);
avgLFP_elD = avgLFP(:,25:48);

%% Get CSD
CSD_elC = calcCSD(avgLFP_elC).*0.4; 
CSD_elD = calcCSD(avgLFP_elD).*0.4;

bl =1:50;
bl_CSD_elC = mean(CSD_elC(:,bl),2);
bl_CSD_elD = mean(CSD_elD(:,bl),2);
blCSD_elC = (CSD_elC-bl_CSD_elC); 
blCSD_elD = (CSD_elD-bl_CSD_elD);

gauss_sigma = 0.1;
padCSD_elC = padarray(blCSD_elC,[1 0],NaN,'replicate');
padCSD_elD = padarray(blCSD_elD,[1 0],NaN,'replicate');
fCSD_elC = filterCSD(padCSD_elC,gauss_sigma);
fCSD_elD = filterCSD(padCSD_elD,gauss_sigma);


%% save fCSD
cd 'E:\LaCie\DATA_KD\BRFS_fCSD'

dCOS_2_elC_161007 = fCSD_elC;
dCOS_2_elD_161007 = fCSD_elD;

saveName = strcat('BRFS_fCSD_dCOS_2_161007');
save(saveName,'dCOS_2_elC_161007','dCOS_2_elD_161007');




