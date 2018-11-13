 function [Cond Unq_cond] = allocateConditions(Grating, pEvC)

a = 0;
b = 0;
c = 0;
d = 0;

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
            if strcmp('Binocular',Grating.stim(t)) %Binocular
                a = a+1;
                    Cond.A(a,1) = Grating.s1_eye(t);
                    Cond.A(a,2) = Grating.s2_eye(t);
                    Cond.A(a,3) = Grating.s1_tilt(t);
                    Cond.A(a,4) = Grating.s2_tilt(t);
                    Cond.A(a,5) = Grating.s1_contrast(t);
                    Cond.A(a,6) = Grating.s2_contrast(t);
                    Cond.A(a,7) = Grating.soa(t);

                    if strcmp('Monocular',Grating.stim(t))
                        Cond.A(a,8) = 1;
                    elseif strcmp('Binocular',Grating.stim(t))
                        Cond.A(a,8) = 2;
                    elseif strcmp('dCOS',Grating.stim(t))
                        Cond.A(a,8) = 3;
                    else
                        Cond.A(a,8) = NaN;
                        disp('error, check "a" for-loop for grating.stim')
                    end

            elseif strcmp('dCOS',Grating.stim(t))
                b = b+1;
                    Cond.B(b,1) = Grating.s1_eye(t);
                    Cond.B(b,2) = Grating.s2_eye(t);
                    Cond.B(b,3) = Grating.s1_tilt(t);
                    Cond.B(b,4) = Grating.s2_tilt(t);
                    Cond.B(b,5) = Grating.s1_contrast(t);
                    Cond.B(b,6) = Grating.s2_contrast(t);
                    Cond.B(b,7) = Grating.soa(t);

                    if strcmp('Monocular',Grating.stim(t))
                        Cond.B(b,8) = 1;
                    elseif strcmp('Binocular',Grating.stim(t))
                        Cond.B(b,8) = 2;
                    elseif strcmp('dCOS',Grating.stim(t))
                        Cond.B(b,8) = 3;
                    else
                        Cond.B(b,8) = NaN;
                        disp('error, check "a" for-loop for grating.stim')
                    end
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
                    Cond.C(c,1) = Grating.s1_eye(t);
                    Cond.C(c,2) = Grating.s2_eye(t);
                    Cond.C(c,3) = Grating.s1_tilt(t);
                    Cond.C(c,4) = Grating.s2_tilt(t);
                    Cond.C(c,5) = Grating.s1_contrast(t);
                    Cond.C(c,6) = Grating.s2_contrast(t);
                    Cond.C(c,7) = Grating.soa(t);

                    if strcmp('Monocular',Grating.stim(t))
                        Cond.C(c,8) = 1;
                    elseif strcmp('Binocular',Grating.stim(t))
                        Cond.C(c,8) = 2;
                    elseif strcmp('dCOS',Grating.stim(t))
                        Cond.C(c,8) = 3;
                    else
                        Cond.C(c,8) = NaN;
                        disp('error, check "a" for-loop for grating.stim')
                    end

            elseif strcmp('dCOS',Grating.stim(t))
                d = d+1;
                    Cond.D(d,1) = Grating.s1_eye(t);
                    Cond.D(d,2) = Grating.s2_eye(t);
                    Cond.D(d,3) = Grating.s1_tilt(t);
                    Cond.D(d,4) = Grating.s2_tilt(t);
                    Cond.D(d,5) = Grating.s1_contrast(t);
                    Cond.D(d,6) = Grating.s2_contrast(t);
                    Cond.D(d,7) = Grating.soa(t);

                    if strcmp('Monocular',Grating.stim(t))
                        Cond.D(d,8) = 1;
                    elseif strcmp('Binocular',Grating.stim(t))
                        Cond.D(d,8) = 2;
                    elseif strcmp('dCOS',Grating.stim(t))
                        Cond.D(d,8) = 3;
                    else
                        Cond.D(d,8) = NaN;
                        disp('error, check "a" for-loop for grating.stim')
                    end

            else
               disp('error, please check EV.C&D loop') 
            end
        else
            disp('error, please check idx loop')
        end
 end


% create cond, a stim presentation x variable types variable that records
% what was displayed on each and every trial.
for e = 1:length(Grating.stim)
    Cond.all(e,1) = Grating.s1_eye(e);
    Cond.all(e,2) = Grating.s2_eye(e);
    Cond.all(e,3) = Grating.s1_tilt(e);
    Cond.all(e,4) = Grating.s2_tilt(e);
    Cond.all(e,5) = Grating.s1_contrast(e);
    Cond.all(e,6) = Grating.s2_contrast(e);
    Cond.all(e,7) = Grating.soa(e);
end
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