function wtnew = basalupdate(basal, tpnno, ts, stdplength)
%basalupdate update weights on basal datastructure
%   for each basal synapse, check whether it has receives a spike recently
% update depends on deltaT (= ts - last spike time). Can use
% basal(tpnno).basalinputs and ts to find deltaT
wtnew = basal(tpnno).basalsynapseweights ; % initialise
for synapseno = 1:basal(tpnno).n_basalinputs
    if basal(tpnno).basalspikeno > 1 % there has been at least one spike
        deltaT = ts - basal(tpnno).basalinputs(basal(tpnno).basalspikeno - 1) ; % deltaT in timesteps
        if (deltaT < stdplength) % last spike within STDP time: NB Strictly less
            % find location in alphafunction
            alphaval = basal(tpnno).alpha_basal(stdplength-deltaT) ;
            wtnew(synapseno) = basal(tpnno).basalsynapseweights(synapseno) + basal(tpnno).wtchange * alphaval ;
            % do we want to change the delay as well?
        else
            % what to do with synapse when it didn't contribute? 
            wtnew(synapseno) = basal(tpnno).basalsynapseweights(synapseno) * 0.98 ;
            % do we want to change the delay as well?
        end
    end
    %  wtnew = 1.1 * basal(tpnno).basalsynapseweights  ;
end % synapse loop
end