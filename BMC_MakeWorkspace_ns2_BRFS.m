%BMC getStarted 180627 MakeWorkspace
%Creates and saves a matlab workspace with variables used in MUA & LFP
%processing and plotting as well as CSD processing and plotting.


clear
close all
addpath('E:\LaCie\MATLAB\helper functions\MLAnalysisOnline\')
%    cd 'E:\LaCie\DATA_KD\161007_E'
%    Filename ='161007_E_brfs001';
   
% % %   cd 'D:\LaCie\DATA_KD\brfs\Brock\151231_E'
% % %   Filename ='151231_E_brfs001';

  cd 'D:\LaCie\DATA_KD\161005_E'
  Filename = '161005_E_brfs001';

ext    = '.gBrfsGratings';
grating = readBRFS([brdrname BRdatafile ext]);

% 1. READ NEV FILE & EXTRACT EVENT CODES/TIMES
NEV = openNEV(strcat(Filename,'.','nev'),'noread','nomat','nosave');
EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes = double(NEV.Data.SerialDigitalIO.TimeStamp); 
%EventTimes = NEV.Data.SerialDigitalIO.TimeStampSec * 1000; 
[pEvC,pEvT] = parsEventCodesML(EventCodes,EventTimes);
evtFs = double(NEV.MetaTags.SampleRes);
clear NEV 
 
% 2. FIND STIMULUS EVENTS/TIMES
obs  = 0;
pre  = 256/1000; % 256ms
post = 612/1000; % 612ms
trls = find(cellfun(@(x) sum(x == 23) == sum(x == 24),pEvC));
for tr = 1:length(trls)
        t = trls(tr);
        
        if ~any(pEvC{t} == 96) % This is not necessary on the evp trails
            % skip if trial was aborted and animal was not rewarded (event code 96)
            continue
        end
 
        stimon1   =  pEvC{t} == 23 | pEvC{t} == 25;
        stimoff1  =  pEvC{t} == 24 | pEvC{t} == 26;
         stimon1   =  pEvC{t} == 23 | pEvC{t} == 25;
        stimoff1  =  pEvC{t} == 24 | pEvC{t} == 26;
        
        start   =  pEvT{t}(stimon1);
        stop    =  pEvT{t}(stimoff1);
        
        maxpres = length(stop);
 
    for p = 1:maxpres
            
            % trigger point
            obs = obs + 1;
            
            % file info
            EV.ec(obs,:)     = [stimon1(p) stimoff1(p)];
            EV.tp(obs,:)     = [start(p) stop(p)];
    end
end
clearvars -except EV Filename
%%
% contrast, soa, orientation, eye

unq_contrast_s1    = nanunique(grating.s1_contrast); 
unq_contrast_s2    = nanunique(grating.s2_contrast); 
unq_ori_s1         = nanunique(grating.s1_tilt); 
unq_ori_s2         = nanunique(grating.s2_tilt); 
unq_soa            = nanunique(grating.soa); 
unq_eye_s1         = nanunique(grating.s1_eye); 
unq_eye_s2         = nanunique(grating.s2_eye);
unq_stim           = unique(grating.stim);
% vec = [1 1; 1 1; 2 2]; 
% unique(vec,'rows')

cond = zeros(length(grating.stim),7 );

for c = 1:length(grating.stim)
    
    cond(c,1) = grating.s1_eye(c);
    cond(c,2) = grating.s2_eye(c);
    cond(c,3) = grating.s1_tilt(c);
    cond(c,4) = grating.s2_tilt(c);
    cond(c,5) = grating.s1_contrast(c);
    cond(c,6) = grating.s2_contrast(c);
    cond(c,7) = grating.soa(c);
    
end
unq_cond = nanunique(cond,'rows');
cellcond = num2cell(cond);

%%% TEST.
%  for i:length(grating.stim) 
%  if grating.stim = 'dCOS'
%       dCOS(i,:) = cond(i,:)
%  else 
%       dCOS(i,:) = nan
%
%

%%
% 3. LOAD MATCHING NEURAL DATA
if exist(strcat(Filename,'.ns2'),'file') == 2
    extension = 'ns2';
end
% 3.1 Read in NS Header
NS_Header = openNSx(strcat(Filename,'.',extension),'noread');

% 3.2 get basic info about recorded data
neural = ~strcmp('E',{NS_Header.ElectrodesInfo.ConnectorBank}); % bank E is the BNCs on the front of the NSP
N.electrodes = length(neural);
N.neural = sum( neural);
N.analog = sum(~neural);

% 3.3 get labels
NeuralLabels = {NS_Header.ElectrodesInfo(neural).Label};
NeuralInfo = NS_Header.ElectrodesInfo(neural);
 
% 3.4 get sampling frequency
Fs = NS_Header.MetaTags.SamplingFreq;

% 3.5 process data electrode by electrode
clear act nct
nct = 0; % prepare counters
for e = 1:N.electrodes
  e
  
    if neural(e) == 1
        nct = nct+1;
        
        clear NS DAT
        electrode = sprintf('c:%u',e);       
        NS = openNSx(strcat(Filename,'.',extension),electrode,'read','uV','precision','double');
        if iscell(NS.Data)
            DAT = cell2mat(NS.Data); 
        else
            DAT = NS.Data;
        end
        NS.Data = [];
                        
        % f_calcMUA
% % % %         mua = f_calcMUA(DAT,Fs,'extralp');
       
        %preallocation
        if nct == 1
% % % %             N.samples.mua = length(mua); 
            N.samples.DAT = length(DAT);
            clear MUA rawMUA LFP rawLFP
% % % %             fMUA = zeros(N.samples.mua,N.neural);
            rawLFP = zeros(ceil(N.samples.DAT),N.neural);
        end
% % % %         fMUA(:,nct) = mua;
        rawLFP(:,nct) = DAT;
        clear mua DAT
        
    end
    
end
clearvars -except fMUA rawLFP EV Filename NeuralLabels Fs extension
%% sort electrodes order
electrodes = unique(NeuralLabels);
 for ch = 1:length(electrodes)
    chname = electrodes{ch}; 
    id = find(~cellfun('isempty',strfind(NeuralLabels,chname)));
    if ~isempty(id)
        ids(ch) = id;
    end
 end  

 clear MUA prefilteredLFP

% % % % MUA = fMUA(:,ids);
clear fMUA
prefilteredLFP = rawLFP(:,ids);
clear rawLFP

%% Filter LFP (MUA was already filtered in 3.5's loop with f_calcMUA function)
% Eventually, this section should be consolidated to a function
lpc = 300; %low pass cutoff
nyq = Fs/2;
lWn = lpc/nyq;
[bwb,bwa] = butter(4,lWn,'low');
for chan = 1:size(prefilteredLFP,2)
    lpLFP(:,chan) = filtfilt(bwb,bwa,prefilteredLFP(:,chan));  %low pass filter to avoid aliasing
end
clear prefilteredLFP bwb bwa 1Wn
clear LFP
LFP = lpLFP; % placeholder. only needed for ns2 version. see below.
% LFP = downsample(lpLFP,30); % downsample to 1kHz %the ns2 should be at
% 1kHz. This step is not necessary in the ns2 version.
clear lpLFP

clearvars -except LFP MUA EV Filename extension

% %Final variable exports are 'MUA' and 'onekHzLFP'.
% %Units are 625532x48 double. This is samples x chan.
% save WorkspaceForProcessing_180705_evp001_test.mat
% 
% %CONGRATS... it finished
% load gong
% sound(y,Fs)

%% save the workspace

%Final variable exports are 'MUA' and LFP'.
%Units are 625532x48 double. This is samples x chan.
dateFormatOut = 'yyyy-mm-dd';
saveDate = datestr(now,dateFormatOut);
saveName = strcat('TESTWorkspace','_',Filename,extension,'_',saveDate);
clear dateFormatOut saveDate
save(saveName);


%notification sound
load gong
sound(y,Fs)



%% 6. HELPER FUNCTIONS
function [pEvC, pEvT] = parsEventCodesML(EventCodes,EventTimes)
 
if isempty(EventTimes) || isempty(EventCodes)
    pEvC = {};
    pEvT = {};
    return
end
 
% these first 2 if statements allow function to be run as:
% parsEventCodesML(NEV.Data.SerialDigitalIO.UnparsedData,NEV.Data.SerialDigitalIO.TimeStampSec)
if ~any(EventCodes == 9)
   EventCodes = EventCodes - 128;
end
if  EventTimes(1)<1
    EventTimes = EventTimes * 1000; %ms, to match 1kHz
end
 
% get trial start index
stind = find(EventCodes == 9);
d = diff(diff(stind)); d = [NaN; d; NaN];
stind = stind(d ==0);
ntr   = length(stind);
if EventCodes(end) ~= 18 % data collection was stoped within a trial
    ntr = ntr-1;
end
 
% setup output vars
pEvC = cell(1,ntr);
pEvT = cell(1,ntr);
 
 
for tr = 1:ntr
    ind = stind(tr)-1;
    ct = 0;
    while EventCodes(ind) ~= 18
        ct = ct + 1;
       
        pEvC{tr}(ct,1) = EventCodes(ind);
        pEvT{tr}(ct,1) = EventTimes(ind);
       
        ind = ind+1;
    end
    % get final 18s
    for ending = 1:3
        pEvC{tr}(ct+ending,1) = EventCodes(ind+ending-1);
        pEvT{tr}(ct+ending,1) = EventTimes(ind+ending-1);
    end
 
end
% convert to double
pEvC = cellfun(@double,pEvC,'UniformOutput',0);
pEvT = cellfun(@double,pEvT,'UniformOutput',0);
end

function MUA = f_calcMUA(DAT,Fs,method)
 
 
if nargin < 3
    method = 'default';
end
 
% Nyquist frequency
nyq = Fs/2;
 
 
switch method
    
    case 'default'
        hpc = 750;  %high pass cutoff
        hWn = hpc/nyq;
        [bwb,bwa] = butter(4,hWn,'high');
        hpMUA = abs(filtfilt(bwb,bwa,DAT)); %high pass filter &rectify
        
        lpc = 200; %low pass cutoff
        lWn = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        lpMUA = filtfilt(bwb,bwa,hpMUA);  %low pass filter to smooth
        
    case 'extralp'
         hpc = 750;  %high pass cutoff
        hWn = hpc/nyq;
        [bwb,bwa] = butter(4,hWn,'high');
        hpMUA = abs(filtfilt(bwb,bwa,DAT)); %high pass filter &rectify
        
        lpc = 50; %low pass cutoff
        lWn = lpc/nyq;
        [bwb,bwa] = butter(4,lWn,'low');
        lpMUA = filtfilt(bwb,bwa,hpMUA);  %low pass filter to smooth
                
end
 
if Fs > 1000
    r = Fs/1000; % 1000 is the sampeling frequency we want after decimation
    MUA(:,1) = downsample(lpMUA,r); % downsample to 1kHz, lowpass filtered previously to avoid aliasing
else
    MUA(:,1) = lpMUA;
end
end


