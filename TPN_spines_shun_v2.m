function [apicalactivation, basalactivation,ahactiv,  threshold, neuron , spikelist] = TPN_spines_shun_v2(varargin)
% TPN_spines_shunt: A version of a two point neuron based on TPN.m,but with spines
% at each synapse
% Dendrite has spines, each with a spine resistance and capacitance and
% leakiness
% Dendrite spines are attached to has its own R, C and leakiness as well
% Introduces shunting synapses, applid to apical and basal activation
%
% All parameters are  <name,value> pairs using varargin
%% (in structure simulation)
% duration: time in seconds for running the neuron
% timestep:  timestep length in seconds
%
%% (in structure apical)
% noapicalinputs:  number of apical inputs
% apicalinputs: actual apical spike inputs in <time synapse_number> format
% apicalsynapseweights: weights on apical synapses, must be same length as noapicalinputs
% tau_apical: time constant for apical synapses (alpha function)
% concept is that dendrite has a local capacitance and membrane
% resistance. In addition there is a resistance describing the maximal
% conductance (=minimum resistance) of the dendrite that governs the amount
% of charge transferred. Keep one set of values for apical and one for
% basal.
% c_apical: capacitance at apical dendrite
% r_apical: apical leakage resistance
% r_synap_dendrite:
% added for spines: currently same for all spines, but this is likely to
% be what gets modified for adaptation
% c_apical_spine: capacitance at apical spinendrite: resistance (1/conductance) of apical dendrite.
% r_apical spine:
% r_synap_spine: resistance (1/conductance) of apical spine. weight
% (apicalsynapseweights) adjusts this per synapse% r_apical_spine: apical spine leakage resistance

%% (in structure basal)
% nobasalinputs: number of basal inputs
% basalinputs: same format
% basalsynapsewights: weights on basal synapses, must be same length as nobasalinputs
% tau-basal: time constant for basal synapses (alpha function)
% concept is that dendrite has a local capacitance and membrane
% resistance. In addition there is a resistance describing the maximal
% conductance (=minimum resistance) of the dendrite that governs the amount
% of charge transferred. Keep one set of values for apical and one for
% basal.
% r_synba_dendrite: resistance of basal dendrite:
% c_basal: capacitance at basal dendrite
% r_basal: basal leakage resistance
% added for spines: currently same for all spines, but this is likely to
% be what gets modified for adaptation
% c_basal_spine: capacitance at basal spine
% r_basal_spine: basal lspineeakage resistance
% r_synba_spine: resistance of basal spine: weight
% (basalsynapseweights) adjusts this per synapse

%% (in structure shunts)
% noapicalshunts: number of apical shunting synapses
% nobasalshunts: number of basal shunting synapses
% apicalshuntinputs: actual apical shunt inputs, in <time synapse_number> format
% basalshuntinputs: actual basal shunt inputs, in <time synapse_number> format
% apicalshuntweights: weights for apical shunt: 0 leq value leq 1
% basalshuntweights: weights for basal shunt: 0 leq value leq 1
% apicalshuntduration: time that apical shunt is active for (in seconds)
% basalshuntduration: time that basal shunt is active for (in seconds)
%
%

%% (in structure neuron)
% parameters for spike generation
% thresh_value: initial threshold for dtetermining spike (applied to ahactiv
% value)
% refractoryperiod: refractory period during which no spikes will be
% generated.
% relrefperiod: relative refractory period after which threshold returns to
%  initial value
% thresh_leap : leap in threshold value after Refractory period
% thresh_decay : exponent of decay (exp(-thresh_decay * time))
% parameters for returning spike list
% neuronid: identity of this neuron (integer)
% maxnospikes: maximum number of spikes this neuron can return
%
% returns apicalactivation, basalactivation
%  axon hillock activity (ahactiv) level and spike times array <neuronid time>
%
% TPN.m started 16 July 2024 by LSS, this new function is based on
% the version 12 8 24 LSS
% TPN_spines started 16 August, LSS.
% TPN_spines_shun started 19 August 2024 LSS
% added apical and basal shunting synapses (working) 3 Sept 2024
% new version allowing loading and saving of , as well as putting the 
% runstep part as an external function. started 24 Sept
% 2024. Parameters now in structures simulation, neuron, apical, basal,
% shunts. Appears working 1 Oct 2024. LSS

% get the parameters using  varargin name value system
%
% default neuronid
neuron.neuronid = 1 ;
% default max number of spikes
neuron.maxnospikes = 1000 ;
% default timestep
simulation.timestep= 0.0001; % 100 microseconds
% lest there be no shunting information, ensure program doesn't crash
shunts.noapicalshunts = 0 ;
shunts.nobasalshunts = 0 ;

%
i=1 ;
while(i<=size(varargin,2))
    switch lower(varargin{i})
        case 'duration' % simulation
            simulation.duration=varargin{i+1};
            i=i+1;
        case 'timestep' % simulation
            simulation.timestep=varargin{i+1};
            i=i+1;
        case 'noapicalinputs' % TPN no of apical inputs
            apical.noapicalinputs =varargin{i+1};
            i=i+1;
        case 'nobasalinputs' % TPN no of basal inputs
            basal.nobasalinputs=varargin{i+1};
            i=i+1;
        case 'apicalinputs' % (this simulation only) actual apical spike inputs
            apical.apicalinputs = varargin{i+1};
            i=i+1;
        case 'basalinputs' % (this simulation only) actual basal spike inputs
            basal.basalinputs = varargin{i+1};
            i=i+1;
        case 'apicalsynapseweights' % TPN weights on apical synapses
            % Must be noapicalinputs long
            apical.apicalsynapseweights = varargin{i+1};
            if (length(apical.apicalsynapseweights) ~= apical.noapicalinputs)
                error("Apical synapse weight vector not same as number of apical inputs") ;
            end
            i=i+1;
        case 'basalsynapsewights' % TPN weights on basal synapses
            % Must be nobasalinputs long
            basal.basalsynapseweights = varargin{i+1};
            if (length(basal.basalsynapseweights) ~= basal.nobasalinputs)
                error("Basal synapse weight vector not same as number of basal inputs") ;
            end
            i=i+1;
        case 'tau_apical' % all TPNs used in alpha function for apical synapses
            apical.tau_apical = varargin{i+1};
            i=i+1;
        case 'tau_basal' % all TPNs used in alpha function for basal synapses
            basal.tau_basal = varargin{i+1};
            i=i+1;
        case 'c_apical' %  all TPNs capacitance of apical dendrite (lumped: used after summing currents)
            apical.c_apical = varargin{i+1};
            i=i+1;
        case 'c_basal' %  all TPNs capacitance of basal dendrite (lumped: used after summing currents)
            basal.c_basal = varargin{i+1};
            i=i+1;
        case 'r_apical' % all TPNs resistance across c_apical: inverse of leak (lumped: used after summing currents)
            apical.r_apical = varargin{i+1};
            i=i+1;
        case 'r_basal' %  all TPNs resisitance across c_basal: inverse of leak (lumped: used after summing currents)
            basal.r_basal = varargin{i+1};
            i=i+1;
        case 'r_synap_dendrite' % all TPNs  resistance of apical synapses
            apical.r_synap_dendrite = varargin{i+1};
            i=i+1;
        case 'r_synba_dendrite' %  all TPNs resistance of basal synapses
            basal.r_synba_dendrite = varargin{i+1};
            i=i+1;
        case 'c_apical_spine' % all TPNs  capacitance of apical dendrite (lumped: used after summing currents)
            apical.c_apical_spine = varargin{i+1};
            i=i+1;
        case 'c_basal_spine' %  all TPNs capacitance of basal dendrite (lumped: used after summing currents)
            basal.c_basal_spine = varargin{i+1};
            i=i+1;
        case 'r_apical_spine' % all TPNs  resistance across c_apical: inverse of leak (lumped: used after summing currents)
            apical.r_apical_spine = varargin{i+1};
            i=i+1;
        case 'r_basal_spine' % all TPNs  resisitance across c_basal: inverse of leak (lumped: used after summing currents)
            basal.r_basal_spine = varargin{i+1};
            i=i+1;
        case 'r_synap_spine' %  all TPNs resistance of apical synapses
            apical.r_synap_spine = varargin{i+1};
            i=i+1;
        case 'r_synba_spine' % all TPNs  resistance of basal synapses
            basal.r_synba_spine = varargin{i+1};
            i=i+1;
        case 'thresh_value' %  all TPNs initial threshold for spiking
            neuron.thresh_value = varargin{i+1};
            i=i+1;
        case 'thresh_leap' %  all TPNs amount by which threshold jumps
            neuron.thresh_leap = varargin{i+1};
            i=i+1;
        case 'thresh_decay' %  all TPNs decay rate for threshold
            neuron.thresh_decay = varargin{i+1};
            i=i+1;
        case 'refractoryperiod' %  all TPNs period during which no spikes can be generated after a spike
            neuron.refractoryperiod = varargin{i+1};
            i=i+1;
        case 'relrefperiod' %  all TPNs relative refractory period
            neuron.relrefperiod = varargin{i+1};
            i=i+1;
        case 'neuronid' % Neuron  identity of neuron: integer
            neuron.neuronid = varargin{i+1};
            i=i+1;
        case 'maxnospikes' % Neuron max number of spikes for this neuron
            neuron.maxnospikes = varargin{i+1};
            i=i+1;
        case 'noapicalshunts' % TPN number of apical shunting synapses
            shunts.noapicalshunts = varargin{i+1};
            i=i+1;
        case 'nobasalshunts' % TPN number of basal shunting synapses
            shunts.nobasalshunts = varargin{i+1};
            i=i+1;
        case 'apicalshuntinputs' % (this simulation only) actual apical shunt inputs, in <time synapse_number> format'
            shunts.apicalshuntinputs = varargin{i+1};
            i=i+1;
        case 'basalshuntinputs' % (this simulation only) actual basal shunt inputs, in <time synapse_number> format'
            shunts.basalshuntinputs = varargin{i+1};
            i=i+1;
        case 'apicalshuntweights' % TPN weights for apical shunt: 0 leq value leq 1'
            shunts.apicalshuntweights = varargin{i+1};
            if (length(shunts.apicalshuntweights) ~= shunts.noapicalshunts)
                error("Apical shunting synapse weight vector not same as number of apical shunts") ;
            end
            i=i+1;
        case 'basalshuntweights' % TPN weights for basal shunt: 0 leq value leq 1'
            shunts.basalshuntweights = varargin{i+1};
            if (length(shunts.basalshuntweights) ~= shunts.nobasalshunts)
                error("Basal shunting synapse weight vector not same as number of basal shunts") ;
            end
            i=i+1;
        case 'apicalshuntduration' % all TPNs  ime that apical shunt is active for (in seconds)'
            shunts.apicalshuntduration = varargin{i+1};
            i=i+1;
        case 'basalshuntduration' % All TPNs time that basal shunt is active for (in seconds)
            shunts.basalshuntduration = varargin{i+1};
            i=i+1;
        otherwise
            error('TPN.m: Unknown argument %s given',varargin{i});
    end
    i=i+1;
end

% length in timesteps of simulation
simulation.simlength = ceil(simulation.duration/simulation.timestep) ;

% calculate alpha functions for apical and basal dendrites (tau need not be
% the same: we might later want a variety of these for excitatory and inhibitory synapses)
apical.alpha_apical = setupAlphaFunction(simulation.timestep, apical.tau_apical) ;
basal.alpha_basal = setupAlphaFunction(simulation.timestep, basal.tau_basal );

% array for apical synapse activations (extra length to stop late spikes
% falling off the end)
% the per-synapse post-tynaptic values are used in this version
apicalpostsynapse = zeros([apical.noapicalinputs simulation.simlength  + length(apical.alpha_apical)]);
% array for basal activations
basalpostsynapse = zeros([basal.nobasalinputs simulation.simlength + length(basal.alpha_basal)]);
% vector for apical current
apicalcurrent = zeros([1 simulation.simlength + length(apical.alpha_apical)]) ;
% vector for basal current
basalcurrent = zeros([1 simulation.simlength + length(basal.alpha_basal)]) ;
% vector for apical activation (charged by apicalcurrent, capacitor and
% parallel resistor)
apicalactivation = zeros([1 simulation.simlength + length(apical.alpha_apical)]) ;
% vector for basal activation (charged by basalcurrent, capacitor and
% parallel resistor)f
basalactivation = zeros([1 simulation.simlength + length(basal.alpha_basal)]) ;
% vector for axon hillock activation
ahactiv = zeros([1 simulation.simlength]) ;
% vector for threshold: useful for seeing what has been going on
% threshold = ones([1 simlength]) * thresh_value ; NO: cannot declare till
% calc-thresh_increment has been calculated
% vector for recording spikes
% spikevector = zeros([1 simlength]) ; no longer required


% sort apical and basal inputs, additive and shunting, into time order
% both apical and basal inputs are a 2d array with
% N_spikes rows 2 columns <time neuron>
if  (isfield(apical, 'apicalinputs') && ~isempty(apical.apicalinputs) )% allow empty apical input
    apical.apicalinputs = sortrows(apical.apicalinputs, 1) ;
    % replace times with timestep numbers
    apical.apicalinputs(:,1) = round(apical.apicalinputs(:,1)/simulation.timestep);
    % check that the apical spike inputs are within range
    if (max(apical.apicalinputs(:,2)) > apical.noapicalinputs)
        error("TPN.m: apical input out of range");
    end
end
if  (isfield(basal, 'basalinputs') && ~isempty(basal.basalinputs)) % allow empty basal input
    basal.basalinputs = sortrows(basal.basalinputs, 1) ;
    % replace times with timestep numbers
    basal.basalinputs(:,1) = round(basal.basalinputs(:,1)/simulation.timestep);
    % check that the basal spike inputs are within range
    if (max(basal.basalinputs(:,2)) > basal.nobasalinputs)
        error("TPN.m: basal input out of range");
    end
end
if (isfield(shunts,'apicalshuntinputs') && ~isempty(shunts.apicalshuntinputs))
    shunts.apicalshuntinputs = sortrows(shunts.apicalshuntinputs, 1);
    % replace times with timestep numbers
    shunts.apicalshuntinputs(:,1) = round(shunts.apicalshuntinputs(:,1)/simulation.timestep);
    % check that the apical spike inputs are within range
    if (max(shunts.apicalshuntinputs(:,2)) > shunts.noapicalshunts)
        error("TPN.m: apical shunt input out of range");
    end
end
if (isfield(shunts, 'basalshuntinputs') && ~isempty(shunts.basalshuntinputs))
    shunts.basalshuntinputs = sortrows(shunts.basalshuntinputs, 1);
    % replace times with timestep numbers
    shunts.basalshuntinputs(:,1) = round(shunts.basalshuntinputs(:,1)/simulation.timestep);
    % check that the basal spike inputs are within range
    if (max(shunts.basalshuntinputs(:,2)) > shunts.nobasalshunts)
        error("TPN.m: basal shunt input out of range");
    end
end

% set up shunting synapses but only if there are shunting synapses
if (shunts.noapicalshunts > 0)
    shunts.asduration = floor(shunts.apicalshuntduration/simulation.timestep) ; % get apical shunt duration in timesteps
    if (shunts.asduration < 1)
        shunts.asduration = 1 ;
    end % 0 is possible, but needs reset to 1
    % calculate shunt value per timestep for each shunting synapse weight
    shunts.ap_sh_wt_pt = zeros([1 length(shunts.apicalshuntweights)]) ; % preallocate
    for sno = 1: length(shunts.apicalshuntweights)
        shunts.ap_sh_wt_pt(sno) = (1-shunts.apicalshuntweights(sno))^(1/shunts.asduration); % weight to apply per timestep
    end
end
if (shunts.nobasalshunts  > 0)
    shunts.bsduration = floor(shunts.basalshuntduration/simulation.timestep) ; % get basal shunt duration in timesteps
    if (shunts.bsduration < 1)
        shunts.bsduration = 1 ;
    end % 0 is possible, but needs reset to 1
    shunts.ba_sh_wt_pt = zeros([1 length(shunts.basalshuntweights)]) ; % preallocate
    for sno = 1: length(shunts.basalshuntweights)
        shunts.ba_sh_wt_pt(sno) = (1 - shunts.basalshuntweights(sno))^(1/shunts.bsduration); % weight to apply per timestep
    end
end


% calculate apical leakiness per timestep:how much leaks away from
% apicalactivation voltage in 1 timestep
% C * dV/dt = I = V/R => dV/dt = V/(R*C)
% so fraction that leaks away is timestep/(R*C)
apical.apicalfracleak = simulation.timestep / (apical.r_apical * apical.c_apical) ;
basal.basalfracleak = simulation.timestep / (basal.r_basal * basal.c_basal) ;
apical.apicalspinefracleak = simulation.timestep/(apical.r_apical_spine  * apical.c_apical_spine) ;
basal.basalspinefracleak = simulation.timestep/(basal.r_basal_spine  * basal.c_basal_spine) ;


% calculate the amount tio be added to the threshold whne a spike occurs.
neuron.thresh_increment = calc_thresh_increment(neuron.thresh_leap, neuron.thresh_decay, ...
    neuron.refractoryperiod, neuron.relrefperiod, simulation.timestep) ;
neuron.th_inc_length = length(neuron.thresh_increment) ;
% vector for threshold: useful for seeing what has been going on
threshold = ones([1 (simulation.simlength + neuron.th_inc_length)]) * neuron.thresh_value ;



% simulate by running step by step
apical.apicalspikeno = 1;  % where are we in the list of apical spikes
basal.basalspikeno = 1;  % where are we in the list of basal spikes

if (isfield(shunts, 'apicalshuntinputs') && ~isempty(shunts.apicalshuntinputs))
    shunts.apicalshuntspikeno = 1; % where we are in list of apical shunting spikes
    shunts.inapicaltimeinterval = zeros([1 size(shunts.apicalshuntinputs, 1)]) ;
    shunts.countdown_ap = zeros([1 size(shunts.apicalshuntinputs, 1)]) ;
end
if (isfield(shunts, 'basalshuntinputs') && ~isempty(shunts.basalshuntinputs))
    shunts.basalshuntspikeno = 1 ; % where we are in list of basal shunting spikes
    shunts.inbasaltimeinterval = zeros([1 size(shunts.basalshuntinputs, 1)]) ;
    shunts.countdown_bs = zeros([1 size(shunts.basalshuntinputs, 1)]) ;
end
no_of_spikes = 0 ;
spikelist = zeros([neuron.maxnospikes 2]) ;
spikelist(: , 1) = neuron.neuronid ;
for tstep = 1:simulation.simlength
    % has a spike been generated? If so, store it
    % spiked = TPN_runstep(tstep) ; % now an external function
    [spiked, neuron, apicalpostsynapse, basalpostsynapse, apicalcurrent, basalcurrent, ...
    apicalactivation, basalactivation, threshold, ...
   ahactiv,  apical, basal, shunts] = TPN_runstep(tstep, neuron, apical, basal, shunts, ... % parameters
   apicalpostsynapse, basalpostsynapse, apicalcurrent, basalcurrent, apicalactivation, basalactivation, ...
   threshold, ahactiv) ;
    if (spiked)
        no_of_spikes = no_of_spikes + 1 ;
        if (no_of_spikes == neuron.maxnospikes)
            error('TPN.m: max number of spikes for neuron %i exceeded.', neuronid);
        end
        spikelist(no_of_spikes, 2) = tstep * simulation.timestep ;
    end

end

spikelist = spikelist(1:no_of_spikes, :) ; % return only the used part

end