function wtnew = apicalupdate(apical, tpnno, ts)
%apicalupdate update apical weights on spiking
%   for each apical synapse, check whether it has received a spike recently
% update depends on deltaT (= ts - last spike time). Can use
% apical(tpnno).apicalinputs and ts to find deltaT
wtnew = apical(tpnno).apicalsynapseweights ; % initialise
for synapseno = 1:apical(tpnno).n_apicalinputs
    if apical(tpnno).apicalspikeno > 1 % there has been at least one spike
        deltaT = ts - apical(tpnno).apicalinputs(apical(tpnno).apicalspikeno - 1) ; % deltaT in timesteps
        % precise way synapse weight TBD
        % for now leave as is
        wtnew(synapseno) = apical(tpnno).apicalsynapseweights(synapseno) ;
    end
end % synapse loop
  % wtnew = apical(tpnno).apicalsynapseweights * 0.9 ;
end