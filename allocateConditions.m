 function [BRFScond Unq_cond] = allocateConditions(Grating, pEvC)

dCOS = 0;
PS = 0;
NSP = 0;


clear Cond

for tr = 1:length(pEvC)
        t = tr;
        
        if ~any(pEvC{t} == 96) % This is not necessary on the evp trails
            % skip if trial was aborted and animal was not rewarded (event code 96)
            continue
        end

        idx = find(pEvC{t} == 23);
        
        % Assign Ev time points with and without soa.
        if      numel(idx) == 1     % there is no soa.
           % create EV.A&B
            if strcmp('dCOS',Grating.stim(bd)) && Grating.soa(bd) == 0
                dCOS = dCOS+1;
                    Cond.dCOS(dCOS).s1_eye     = Grating.s1_eye(t);
                    Cond.dCOS(dCOS).s2_eye     = Grating.s2_eye(t);
                    Cond.dCOS(dCOS).s1_tilt    = Grating.s1_tilt(t);
                    Cond.dCOS(dCOS).s2_tilt	= Grating.s2_tilt(t);
                    Cond.dCOS(dCOS).s1_contrast = Grating.s1_contrast(t);
                    Cond.dCOS(dCOS).s2_contrast = Grating.s2_contrast(t);
                    Cond.dCOS(dCOS).soa        = Grating.soa(t);
                    Cond.dCOS(dCOS).trial      = Grating.trial(t);
                    Cond.dCOS(dCOS).timestamp = Grating.timestamp(bd);
                    

            elseif strcmp('dCOS',Grating.stim(bd)) && Grating.soa(bd) == 0
                dCOS = dCOS+1;
                    Cond.dCOS(dCOS).s1_eye     = Grating.s1_eye(t);
                    Cond.dCOS(dCOS).s2_eye     = Grating.s2_eye(t);
                    Cond.dCOS(dCOS).s1_tilt    = Grating.s1_tilt(t);
                    Cond.dCOS(dCOS).s2_tilt	= Grating.s2_tilt(t);
                    Cond.dCOS(dCOS).s1_contrast = Grating.s1_contrast(t);
                    Cond.dCOS(dCOS).s2_contrast = Grating.s2_contrast(t);
                    Cond.dCOS(dCOS).soa        = Grating.soa(t);
                    Cond.dCOS(dCOS).trial      = Grating.trial(t);
                    Cond.dCOS(dCOS).timestamp = Grating.timestamp(bd);

                    
            elseif strcmp
% Grating.stim is a cell field containing strings. This is reformatted into
% double format where 1=Mo, 2=Bi, 3=dCOS.
%
for gs = 1:length(Grating.stim)
    if strcmp('Monocular',Grating.stim(gs))
        Cond.all(gs,8) = 1;
    elseif strcmp('Binocular',Grating.stim(gs))
        Cond.all(gs,8) = 2;
    elseif strcmp('dCOS',Grating.stim(gs))
        Cond.all(gs,8) = 3;
    else
        Cond(gs,8) = NaN;
        disp('error, check "gs" for-loop for grating.stim')
    end
end
% Finds all possible combinations of stimuli. Should be 64. 
Unq_cond = nanunique(Cond.all,'rows');

 end