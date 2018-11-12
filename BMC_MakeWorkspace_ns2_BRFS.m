%BMC getStarted 180627 MakeWorkspace
%Creates and saves a matlab workspace with variables used inLFP
%processing and plotting as well as CSD processing and plotting.
%181107 test

clear
close all

if ispc
addpath('E:\LaCie\MATLAB\helper functions\MLAnalysisOnline\')
addpath('E:\LaCie\MATLAB\helper functions\MLAnalysisOnline\NPMK-master\NPMK\')
else
addpath('/Volumes/PassportForMac/MATLAB/functions/helper functions/MLAnalysisOnline/NPMK 2.5.1')
addpath('/Volumes/PassportForMac/MATLAB/functions/helper functions/MLAnalysisOnline');
end

ext    = '.gBrfsGratings';
if ispc
brdrname = 'E:\LaCie\DATA_KD\161005_E\';
else
    brdrname = '/Volumes/PassportForMac/DATA_KD/161005_E/';
end
BRdatafile = '161005_E_brfs001';
  cd(brdrname)
  Filename = BRdatafile;
  
ext    = '.gBrfsGratings';
Grating = readBRFS([brdrname BRdatafile ext]);


% 1. READ NEV FILE & EXTRACT EVENT CODES/TIMES
if ispc
   fname = strcat(Filename,'.','nev');
    NEV = openNEV(fname,'noread','nomat','nosave'); %check NPMK version info
else

NEV = openNEV('noread','nomat','nosave'); %check NPMK version info
end
EventCodes = NEV.Data.SerialDigitalIO.UnparsedData - 128;
EventTimes = double(NEV.Data.SerialDigitalIO.TimeStamp); 
%EventTimes = NEV.Data.SerialDigitalIO.TimeStampSec * 1000; 
[pEvC,pEvT] = parsEventCodesML(EventCodes,EventTimes);
evtFs = double(NEV.MetaTags.SampleRes);
clear NEV 
 
% 2. FIND STIMULUS EVENTS/TIMES
obs  = 0;
a = 0;
b = 0;
c = 0;
d = 0;
pre  = 256/1000; % 256ms
post = 612/1000; % 612ms

for tr = 1:length(pEvC)
        t = tr;
        
        if ~any(pEvC{t} == 96) % This is not necessary on the evp trails
            % skip if trial was aborted and animal was not rewarded (event code 96)
            continue
        end
        
        % Logical index of the pEvC field asigned to trial 't' that alings
        % with either start or stop of a stim.
        stimon   =  pEvC{t} == 23;
        stimoff  =  pEvC{t} == 24;
        
        %Based on the logical stimon index previously created, determine if
        %there is soa or no soa. No soa gets labeled as start_noSoa - this
        %will go into groups A and B. With Soa present, start1 and start2
        %are created. - this will go onti groups C and D. This seperation
        %may not be necessary and should be reviewed later. 
        idx = find(stimon);
        if      numel(idx) == 1     % there is no soa.
            start_noSoa  =  pEvT{t}(stimon);
            
        elseif  numel(idx) == 2     %there is indeed soa
            start1  =  pEvT{t}(stimon(idx(1)));
            start2  =  pEvT{t}(stimon(idx(2)));
        else
            disp('error, please check idx loop')
        end
       
        %The time at which one or both stimuli were removed
        stop    =  pEvT{t}(stimoff);
        
% % % % %         % trigger point
% % % % %         obs = obs + 1; 
        
        % Assign Ev time points with and without soa.
        if      numel(idx) == 1     % there is no soa.
           % create EV.A&B
            if strcmp('Binocular',Grating.stim(t)) %Binocular
                a = a+1;
                EV.A(a,:) = [start_noSoa stop];
            elseif strcmp('dCOS',Grating.stim(t))
                b = b+1;
                EV.B(b,:) = [start_noSoa stop];
            elseif strcmp('Monocular', Grating.stim(t)) % skip Monocular
                continue % skips and dismisses all Monocular trials.
                
            else
               disp('error, please check EV.A&B loop') 
               disp(t)
            end
% % % % % %             EV.tpNoSoa(obs,:)   = [start_noSoa stop];            
        elseif  numel(idx) == 2     %there is indeed soa
            % create EV.C&D 
            if strcmp('Binocular',Grating.stim(t)) %Binocular
                c = c+1;
                EV.C(c,:) = [start1 start2 stop];
            elseif strcmp('dCOS',Grating.stim(t))
                d = d+1;
                EV.D(d,:) = [start1 start2 stop];
            else
               disp('error, please check EV.C&D loop') 
             end
% % % % % %             EV.tpSoa(obs,:)     = [start1 start1 stop];
        else
            disp('error, please check idx loop')
        end
    
 
 end


%% organize stimuli conditions
% Ignore Monocular stim for now. Look at dCOS and Binocular stim under
% simultaneous and stimulus onset asynchrony conditions. 
% 4 groups. 
% A == Bi, soa=0; B == dCOS, soa=0; C == Bi, soa=800; D == dCOS, soa=800.


% to look at all condition parameter options
unq_contrast_s1    = nanunique(Grating.s1_contrast); 
unq_contrast_s2    = nanunique(Grating.s2_contrast); 
unq_ori_s1         = nanunique(Grating.s1_tilt); 
unq_ori_s2         = nanunique(Grating.s2_tilt); 
unq_soa            = nanunique(Grating.soa); 
unq_eye_s1         = nanunique(Grating.s1_eye); 
unq_eye_s2         = nanunique(Grating.s2_eye);
unq_stim           = unique(Grating.stim); %note, this is a cell field

%pre-allocate 
Cond = zeros(length(Grating.stim),8 );
% create cond, a stim presentation x variable types variable that records
% what was displayed on each and every trial.
for c = 1:length(Grating.stim)
    Cond(c,1) = Grating.s1_eye(c);
    Cond(c,2) = Grating.s2_eye(c);
    Cond(c,3) = Grating.s1_tilt(c);
    Cond(c,4) = Grating.s2_tilt(c);
    Cond(c,5) = Grating.s1_contrast(c);
    Cond(c,6) = Grating.s2_contrast(c);
    Cond(c,7) = Grating.soa(c);
end
% Grating.stim is a cell field containing strings. This is reformatted into
% double format where 1=Mo, 2=Bi, 3=dCOS.
%
for gs = 1:length(Grating.stim)
    if strcmp('Monocular',Grating.stim(gs))
        Cond(gs,8) = 1;
    elseif strcmp('Binocular',Grating.stim(gs))
        Cond(gs,8) = 2;
    elseif strcmp('dCOS',Grating.stim(gs))
        Cond(gs,8) = 3;
    else
        Cond(gs,8) = NaN;
        disp('error, check "gs" for-loop for grating.stim')
    end
end
% Finds all possible combinations of stimuli. Should be 64. 
Unq_cond = nanunique(Cond,'rows');

clearvars -except EV Filename Grating Cond Unq_cond
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
clearvars -except fMUA rawLFP EV Filename NeuralLabels Fs extension Cond Unq_cond
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

clearvars -except LFP MUA EV Filename extension Cond Unq_cond

% %Final variable exports are 'MUA' and 'LFP'.
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
saveName = strcat('BRFSWorkspace','_',Filename,extension,'_',saveDate);
clear dateFormatOut saveDate
save(saveName);


%notification sound
load gong
sound(y,Fs)

clear y Fs
