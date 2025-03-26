function wtnew = basalupdate(basal, tpnno, ts)
%basalupdate update weights on basal datastructure
%   for each basal synapse, check whether it has receives a spike recently
% update depends on deltaT (= ts - last spike time). Can use
% basal(tpnno).basalinputs and ts to find deltaT
wtnew = basal(tpnno).basalsynapseweights ; % initialise
for synapseno = 1:basal(tpnno).n_basalinputs
    if basal(tpnno).basalspikeno > 1 % there has been at least one spike
    deltaT = ts - basal(tpnno).basalinputs(basal(tpnno).basalspikeno - 1) ; % deltaT in timesteps
        % precise way synapse weight
        % (basal(tpnno).basalsynapseweights(synapseno) is altered TBD
        % for now leave as is
        wtnew(synapseno) = basal(tpnno).basalsynapseweights(synapseno) ;
    end
   %  wtnew = 1.1 * basal(tpnno).basalsynapseweights  ;
end % synapse loop
end