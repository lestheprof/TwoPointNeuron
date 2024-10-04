function [isspike, neuron,  apicalpostsynapse, basalpostsynapse, apicalcurrent, basalcurrent, ...
    apicalactivation, basalactivation, threshold,...
    ahactiv,  apical, basal, shunts] = TPN_runstep(ts, neuron, apical, basal, shunts, ... % parameters
    apicalpostsynapse, basalpostsynapse, apicalcurrent, basalcurrent, apicalactivation, basalactivation, ...
    threshold, ahactiv)
% run simulation for one timestep, returning true if a spike is
% generated
% calculate apical synaptic depolarisation: fix, there can be > 1
% spike in a single timestep so jse while
% new version: standalone function, instead of being inside an
% enclosing function: started 25 Sept 2024, appears working 1 Oct 2024.
%
%% parameters
% ts: timestep,
% neuron, apical, basal, shunts parameters
% apicalpostsynapse, basalpostsynapse, apicalcurrent, basalcurrent, apicalactivation, basalactivation, 
% threshold, ahactiv: supplied (and return) values for each timestep describing values
% inside neuron.
%% returns
% isspike: true if there is a spike at this timestep
% neuron:some parameters that may be altered are returned
% apicalpostsynapse, basalpostsynapse, apicalcurrent, basalcurrent, 
%   apicalactivation, basalactivation, threshold,...
%   ahactiv:  returned values for variables inside neuron
% apical, basal, shunts: return parameters that may be changed

while ((apical.apicalspikeno <= size(apical.apicalinputs,1)) && ...
        (apical.apicalinputs(apical.apicalspikeno,1) == ts)) % calculate for each synapse: is there a new apical spike input?
    % apical spike detected: add alpha function * weight
    % calculate for each synapse used in this version
    apicalpostsynapse(apical.apicalinputs(apical.apicalspikeno,2), ts:(ts + length(apical.alpha_apical) -1)) = ...
        (apicalpostsynapse(apical.apicalinputs(apical.apicalspikeno,2), ts:(ts + length(apical.alpha_apical) -1)) * (1-apical.apicalspinefracleak)) + ...
        (apical.alpha_apical * apical.apicalsynapseweights(apical.apicalinputs(apical.apicalspikeno,2)) / (apical.r_synap_spine + apical.r_synap_dendrite) );
    % adjust using spine capacitance leakage and resistance

    % add this change to the total depolarisation (currently no
    % geometry on apical dendrite)
    apicalcurrent(ts:(ts + length(apical.alpha_apical) -1)) = ...
        apicalcurrent(ts:(ts + length(apical.alpha_apical) -1)) + ...
        apicalpostsynapse(apical.apicalinputs(apical.apicalspikeno,2), ts:(ts + length(apical.alpha_apical) -1) ) ;
    % (alpha_apical * apicalsynapseweights(apicalinputs(apicalspikeno,2)) / r_synap_dendrite )  ;
    apical.apicalspikeno = apical.apicalspikeno + 1 ;
end

% calculate apical activation:
% total charge is weight coulombs!!! That's why the voltages are
% so high: fixed.
% is the apical current a current (I) or a charge (I * Delta T = Q)?
if (ts > 1)
    apicalactivation(ts) = apicalactivation(ts-1) * (1 - apical.apicalfracleak) + ...
        (apicalcurrent(ts)/apical.c_apical) ;
    % ((timestep * apicalcurrent(ts))/c_apical) ;
    % (timestep * (apicalcurrent(ts)/c_apical)) ;
end
% apply apical shunt, if any.
if ((shunts.noapicalshunts > 0) && (isfield(shunts, 'apicalshuntinputs')))
    %   apply shunting synaptic inhibition here.
    %   for each apical shunting synapse, add up the inhibition and then apply.
    % is there a shunting synapse at this time, or has there been one
    % within the last asduration timesteps?
    while (shunts.apicalshuntspikeno <= size(shunts.apicalshuntinputs, 1) && ...
            (shunts.apicalshuntinputs(shunts.apicalshuntspikeno,1) == ts))
        if (shunts.apicalshuntinputs(shunts.apicalshuntspikeno,1) == ts) % if there's a new shunting input, initialise countdown
            shunts.inapicaltimeinterval(shunts.apicalshuntspikeno)  = true;
            % keep the id of the spike synapse number, apicalshuntinputs(apicalshuntspikeno,2)
            % and use it when indexing ap_sh_wt_pt
            shunts.countdown_ap(shunts.apicalshuntspikeno) = shunts.asduration ;
        end
        shunts.apicalshuntspikeno = shunts.apicalshuntspikeno + 1;
    end % while
    % calculate total shunting effect
    shunteffect = 1 ;
    for apsspikeindex = 1:length(shunts.inapicaltimeinterval)
        if (shunts.inapicaltimeinterval(apsspikeindex)  == true)
            shunteffect = shunteffect * shunts.ap_sh_wt_pt(shunts.apicalshuntinputs(shunts.apicalshuntspikeno - 1,2)) ;
            % shunteffect = shunteffect * ap_sh_wt_pt(apsspikeindex) ; % issue: ap_sh_wt_pt see above.
            shunts.countdown_ap(apsspikeindex) = shunts.countdown_ap(apsspikeindex) - 1;
            if (shunts.countdown_ap(apsspikeindex) == 0)
                shunts.inapicaltimeinterval(apsspikeindex) = false ;
            end % if
        end% if
    end  % for
    % apply shunteffect
    apicalactivation(ts) = apicalactivation(ts) * shunteffect ;
end % if



% calculate basal synaptic depolarisation
while ((basal.basalspikeno <= size(basal.basalinputs,1)) ...
        && (basal.basalinputs(basal.basalspikeno,1) == ts))
    % basal spike detected: add alpha function * weight
    % calculate for ecach synapse (not used in this version)
    basalpostsynapse(basal.basalinputs(basal.basalspikeno,2), ts:(ts + length(basal.alpha_basal) -1)) = ...
        (basalpostsynapse(basal.basalinputs(basal.basalspikeno,2), ts:(ts + length(basal.alpha_basal) -1)) *  (1-basal.basalspinefracleak)) + ...
        (basal.alpha_basal * basal.basalsynapseweights(basal.basalinputs(basal.basalspikeno,2)) / (basal.r_synba_spine + basal.r_synba_dendrite)) ;
    % add this change to the total depolarisation
    basalcurrent(ts:(ts + length(basal.alpha_basal) -1)) = ...
        basalcurrent(ts:(ts + length(basal.alpha_basal) -1)) + ...
        basalpostsynapse(basal.basalinputs(basal.basalspikeno,2), ts:(ts + length(basal.alpha_basal) -1)) ;
    % (alpha_basal * basalsynapseweights(basalinputs(basalspikeno,2))  / r_synba_dendrite) ;
    basal.basalspikeno = basal.basalspikeno + 1 ;
end

% calculate basal activation:
% is the basal current a current (I) or a charge (I * Delta T = Q)?
if (ts > 1)
    basalactivation(ts) = basalactivation(ts-1) * (1 - basal.basalfracleak) + ...
        (basalcurrent(ts)/basal.c_basal) ;
    % ((timestep * basalcurrent(ts))/c_basal) ;
    % (timestep * (basalcurrent(timestep)/c_basal)) ;
end
if ((shunts.nobasalshunts  > 0)  && (isfield(shunts, 'basalshuntinputs')))
    %   apply shunting synaptic inhibition here.
    %   for each basal shunting synapse, add up the inhibition and then apply
    %   apply shunting synaptic inhibition here.
    % is there a shunting synapse at this time, or has there been one
    % within the last asduration timesteps?
    while (shunts.basalshuntspikeno <= size(shunts.basalshuntinputs, 1) && ...
            (shunts.basalshuntinputs(shunts.basalshuntspikeno,1) == ts))
        if (shunts.basalshuntinputs(shunts.basalshuntspikeno,1) == ts) % if there's a new shunting input, initialise countdown
            shunts.inbasaltimeinterval(shunts.basalshuntspikeno)  = true;
            % keep the id of the spike synapse number, basalshuntinputs(apicalshuntspikeno,2)
            % and use it when indexing ap_sh_wt_pt
            shunts.countdown_bs(shunts.basalshuntspikeno) = shunts.bsduration ;
        end
        shunts.basalshuntspikeno = shunts.basalshuntspikeno + 1;
    end % while
    % calculate total shunting effect
    shunteffect = 1 ;
    for bssspikeindex = 1:length(shunts.inbasaltimeinterval)
        if (shunts.inbasaltimeinterval(bssspikeindex)  == true)
            shunteffect = shunteffect * shunts.ba_sh_wt_pt(shunts.basalshuntinputs(shunts.basalshuntspikeno - 1,2)) ;
            shunts.countdown_bs(bssspikeindex) = shunts.countdown_bs(bssspikeindex) - 1;
            if (shunts.countdown_bs(bssspikeindex) == 0)
                shunts.inbasaltimeinterval(bssspikeindex) = false ;
            end % if
        end% if
    end  % for
    % apply shunteffect
    basalactivation(ts) = basalactivation(ts) * shunteffect ;


end

% Try equation from Reza & Ahsan (R^2 + 2*R +2*C *
% (1+mod(R)))for axon hillock activation
ahactiv(ts) = basalactivation(ts)^2  +  2 *  basalactivation(ts) +...
    2 * apicalactivation(ts) * (1 + abs(basalactivation(ts))) ;

% decide whether to spike at this timestep, and store spike if so.
if (ahactiv(ts) > threshold(ts)) % spike!
    isspike = 1 ;
    % update threshold
    threshold(ts:ts + neuron.th_inc_length -1) = threshold(ts:ts + neuron.th_inc_length -1) + ...
        neuron.thresh_increment ;
else
    isspike = 0 ;
end
end % runstep
