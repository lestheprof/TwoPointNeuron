function wtnew = apicalupdate(apical, tpnno, ts, stdplength)
%apicalupdate update apical weights on spiking
%   for each apical synapse, check whether it has received a spike recently
% update depends on deltaT (= ts - last spike time). Can use
% apical(tpnno).apicalinputs and ts to find deltaT
wtnew = apical(tpnno).apicalsynapseweights ; % initialise
for synapseno = 1:apical(tpnno).n_apicalinputs
    if apical(tpnno).apicalspikeno > 1 % there has been at least one spike
        deltaT = ts - apical(tpnno).apicalinputs(apical(tpnno).apicalspikeno - 1) ; % deltaT in timesteps
        if (deltaT < stdplength) % last spike within STDP time: N.B strictly less.
            % find location in alphafunction
            alphaval = apical(tpnno).alpha_apical(stdplength-deltaT) ; 
            % arithmetical (unbounded) increase
            wtnew(synapseno) = apical(tpnno).apicalsynapseweights(synapseno) + apical(tpnno).wtchange * alphaval ;
            % do we want to change the delay as well?
        else
            % what to do with synapse when it didn't contribute?
            % geometrical  decrease towards 0
            wtnew(synapseno) = apical(tpnno).apicalsynapseweights(synapseno) * 0.98 ;
            % do we want to change the delay as well?
        end
    end
end % synapse loop
% wtnew = apical(tpnno).apicalsynapseweights * 0.9 ;
end