function spikelist = runTPNIIsimulation(simfile, apicalfile, basalfile, shuntfile, neuronfile, iifile, weightfile, ...
    connectionfile, externalinputs, drivinginputinfo ...
    )
% runTPNIIsimulation runs the simulation.
%
% parameters are fhe file names for all the parameters
% simfile whole simulation values, including number of TPNs and IIs, duration and timestep
% apicalfile info about apical part of each TPN
% basalfile info about basal part of each TPN
% shuntfile infor about shunts on TPNs
% neuronfile info about each TPN
% iifile info about each inhibitory interneuron
% weightfile info about all the weights
% connectionfile info about all the connections between neurons
% externalinputs info about all the externally generated spikes input.
% inputdir has the directory in which the files to be processed reside
% filelist is the list of files to be processed.

% LSS 18 Dec 2024.
% LSS 15 Jan2025: calls setupweights before setupNetworkV2. Matters because
% setupNetworkV2 uses shunting weights to set up values inside shunts
% structure. Also added resetvalue to apical, basal and IIneuron structures
% to optionally allow resetting activation after a spike is generated
% LSS April 6 2025; add input units (started); Sept 3 2025 continued
% LSS Sept 14 2025; add saving of possibly modified weights to file,
% enabling re-use.

saveparamsandarrays = true ;
spikeplot = true ; % plot spikes
spikesave = true ; % save spikes in form <neuron time>
weightsave = true ; % save (amended) weights
networksave = true ; % save network (connectionfile) (includes delays)

% read simulation parameters in from file
[~, simulation] = readnetwork(simfile) ;

% below put values  into structures. Each one is an array of N_TPNs
% structures.
% read apical file from apicalfile.txt
% apical parameters
apical = readapicalfile(apicalfile, simulation) ;
% basal parameters
basal = readbasalfile(basalfile, simulation) ;
% shunt parameters
shunts = readshuntfile(shuntfile, simulation) ;
% TPN neuron parameters
neuron = readneuronfile(neuronfile, simulation) ;


% these are read from a file
[basalinputs, apicalinputs, apicalshuntinputs, basalshuntinputs, IIinputs, Xinputs] = readexternalinput(externalinputs) ;

% read the driving inputs (test and train eventually) in to a structure
% [traindata, testdata] = readaudioinputspikes(drivinginputinfo) ; till
% later


% initialise rest of structures based on parameters, number of inputs and shunts and external inputs
% allocate last one first: essentially pre-allocating

% input neurons
for inputnno = simulation.N_Inputs:-1:1
    inputneuron(inputnno).number = inputnno ; % really just a placeholder
    if isempty(Xinputs)
        inputneuron(inputnno).inputindex = 0 ; % 0 means no input on this input neuron
    else
        inputneuron(inputnno).xinputs = Xinputs(Xinputs(:,1) == inputnno, :) ;
        inputneuron(inputnno).inputindex = 1 ;
    end
end
for inputnno = simulation.N_Inputs:-1:1
    if isempty(inputneuron(inputnno).xinputs)
        inputneuron(inputnno).inputindex = 0 ;
    end
end
%
%% note that weights have been declared but not set

% need to read in weights before calling setupnetwork

% external II inputs format: <II_number time synapse_number>
IIneuron = readIIneuronfile(iifile, simulation) ;

% read in the weights
% weights are per synapse: <neuron_type neuron_number syn_type
% syn_number weight>
% weightfile = "weights1.txt" ;
[basal, apical, shunts, IIneuron] = setupweights(weightfile, basal, apical, shunts, IIneuron) ;
for tpnno = 1:simulation.N_TPNs
    % calculate the amount tio be added to the threshold whne a spike occurs.
    neuron(tpnno).thresh_increment = calc_thresh_increment(neuron(tpnno).thresh_leap, neuron(tpnno).thresh_decay, ...
        neuron(tpnno).refractoryperiod, neuron(tpnno).relrefperiod, simulation.timestep) ;
end

% now set up interconnection
% connectionfile = "network1.txt" ;
% connectionfile has table for interconnection, format described in setupinterconnection
[inputneuron, neuron, IIneuron] = setupinterconnection( simulation, inputneuron, neuron, IIneuron, connectionfile) ;

% process inputneurons: add input spikes to appropriate TPN apical and
% basal structures. This is practical because the external inputs are
% fixed.
apin = apicalinputs ;
bain = basalinputs ;

aind = size(apin, 1) ;
bind = size(bain, 1) ;

%% redefine apicalinputs and basalinputs to include externalinput from input neurons
% count how many there are
% apical and basal
aindex = aind ;
bindex = bind ;
for xno = 1:simulation.N_Inputs
    for exinputno = 1: size(inputneuron(xno).xinputs, 1 )
        % some of these may fo to more than one place.
        for targetno = 1:length(inputneuron(xno).targets)
            if  char(inputneuron(xno).targets(targetno).to_syntype) == 'A'
                aindex = aindex + 1 ;
            else
                if  char(inputneuron(xno).targets(targetno).to_syntype) == 'B'
                    bindex = bindex + 1 ;
                end
            end
        end
    end
end
apicalinputs = zeros([aindex 3]) ; % neuronno time synapseno
apicalinputs(1:aind, :) = apin ;
basalinputs = zeros([bindex 3]) ; % neuronno time synapseno
basalinputs(1:bind) = bain ;

for xno = 1:simulation.N_Inputs
    for exinputno = 1: size(inputneuron(xno).xinputs, 1 )
        % each input gets appended to the target
        for targetno = 1:length(inputneuron(xno).targets)
            switch char(inputneuron(xno).targets(targetno).to_syntype)
                case 'A' % input to an apical part of a TPN
                    % update
                    aind = aind + 1 ;
                    %  add delay
                    stime = inputneuron(xno).xinputs(exinputno, 2) + inputneuron(xno).targets(targetno).delay ;
                    % inputneuron(xno).xinputs(exinputno, 2) = inputneuron(xno).xinputs(exinputno, 2) /simulation.timestep   ;
                    apicalinputs(aind, 2) =  stime ; %time
                    apicalinputs(aind, 1) = inputneuron(xno).targets(targetno).to_nno ; % neuron no
                    apicalinputs(aind, 3) = inputneuron(xno).targets(targetno).to_synno ; % synapse no

                case 'B' % input to a basal part of a TPN
                    % update
                    bind = bind + 1 ;
                    %  add delay
                    stime = inputneuron(xno).xinputs(exinputno, 2) + inputneuron(xno).targets(targetno).delay ;
                    % inputneuron(xno).xinputs(exinputno, 2) = inputneuron(xno).xinputs(exinputno, 2) /simulation.timestep   ;
                    basalinputs(bind, 2) =  stime ; %time
                    basalinputs(bind, 1) = inputneuron(xno).targets(targetno).to_nno ; % neuron no
                    basalinputs(bind, 3) = inputneuron(xno).targets(targetno).to_synno ; % synapse no
            end
        end

    end
end




for tpnno = simulation.N_TPNs:-1:1 % place in correct structure
    basal(tpnno).basalinputs = basalinputs(basalinputs(:,1) == tpnno, :) ;
    basal(tpnno).basalinputs = basal(tpnno).basalinputs(:, 2:3) ;
    apical(tpnno).apicalinputs = apicalinputs(apicalinputs(:,1) == tpnno, :)  ;
    apical(tpnno).apicalinputs = apical(tpnno).apicalinputs(:, 2:3) ;
    shunts(tpnno).apicalshuntinputs = apicalshuntinputs(apicalshuntinputs(:,1) == tpnno, :)  ;
    shunts(tpnno).apicalshuntinputs = shunts(tpnno).apicalshuntinputs(:,2:3);
    shunts(tpnno).basalshuntinputs = basalshuntinputs(basalshuntinputs(:,1) == tpnno, :)  ;
    shunts(tpnno).basalshuntinputs = shunts(tpnno).basalshuntinputs(:,2:3) ;
    % set up weight vectors (now in file reading function)

    % calculate the amount tio be added to the threshold whne a spike occurs.
    % neuron(tpnno).thresh_increment = calc_thresh_increment(neuron(tpnno).thresh_leap, neuron(tpnno).thresh_decay, ...
    %     neuron(tpnno).refractoryperiod, neuron(tpnno).relrefperiod, simulation.timestep) ;
    neuron(tpnno).th_inc_length = length(neuron(tpnno).thresh_increment) ;
    neuron(tpnno).spikes = zeros([1 neuron(tpnno).maxnospikes]) ;
    neuron(tpnno).spikecount = 0 ;
    apical(tpnno).apicalspikeno = 1;  % where are we in the list of apical spikes
    basal(tpnno).basalspikeno = 1;  % where are we in the list of basal spikes
    % now we have possibly multiple neurons causing spikes in shunting
    % synapses we need to initialise the relevant variables
    %if (isfield(shunts, 'apicalshuntinputs') && ~isempty(shunts(tpnno).apicalshuntinputs))
    shunts(tpnno).apicalshuntspikeno = 1; % where we are in list of apical shunting spikes
    shunts(tpnno).inapicaltimeinterval = zeros([1 size(shunts(tpnno).apicalshuntinputs, 1)]) ;
    shunts(tpnno).countdown_ap = zeros([1 size(shunts(tpnno).apicalshuntinputs, 1)]) ;
    %end
    %if (isfield(shunts, 'basalshuntinputs') && ~isempty(shunts(tpnno).basalshuntinputs))
    shunts(tpnno).basalshuntspikeno = 1 ; % where we are in list of basal shunting spikes
    shunts(tpnno).inbasaltimeinterval = zeros([1 size(shunts(tpnno).basalshuntinputs, 1)]) ;
    shunts(tpnno).countdown_bs = zeros([1 size(shunts(tpnno).basalshuntinputs, 1)]) ;
    %end

end

% II (LIF) neuron setup

for IIno = simulation.N_IIs:-1:1 % allocate last one first: essentially pre-allocating
    IIneuron(IIno).thresh_increment = calc_thresh_increment(IIneuron(IIno).thresh_leap, IIneuron(IIno).thresh_decay, ...
        IIneuron(IIno).refractoryperiod, IIneuron(IIno).relrefperiod, simulation.timestep) ;
    IIneuron(IIno).th_inc_length = length(IIneuron(IIno).thresh_increment) ;
    IIneuron(IIno).spikecount = 0 ;
    % transfer the IIinputs to the appropriate IIneuron
    IIneuron(IIno).inputs = IIinputs(IIinputs(:,1) == IIno, :) ;
    IIneuron(IIno).inputs = IIneuron(IIno).inputs(:, 2:3) ;
    IIneuron(IIno).spikeno = 1 ;
end

[simulation, neuron, basal,apical, shunts, apicalcurrent, basalcurrent, ...
    apicalactivation, basalactivation, ahactiv, IIneuron] = ...
    setupnetworkV2(simulation,neuron,basal, apical,shunts, IIneuron) ;

if (saveparamsandarrays)
    save(strcat("paramsarrays",string(datetime("today")), ".mat"), "simulation", "neuron", "basal","apical", "shunts", ...
        "apicalcurrent", "basalcurrent", ...
        "apicalactivation", "basalactivation", "ahactiv", "IIneuron", "inputneuron");
end


%% run simulation: run TPN_runstep and II_runstep for each TPN and each inhibitory interneuron

for ts = 1:simulation.simlength
    % process TPNs
    for tpnno = 1:simulation.N_TPNs
        [TPNspike, neuron,   apicalcurrent, basalcurrent, apicalactivation, basalactivation, ...
            ahactiv,  apical, basal, shunts] = TPN_runstep(ts, tpnno, simulation, neuron, inputneuron, apical, basal, shunts, ... % parameters
            apicalcurrent, basalcurrent, apicalactivation, basalactivation, ahactiv) ;
        if TPNspike
            if (isfield(neuron(tpnno),"targets"))
                % process spike by supplying spikes to neuron(tpnno).targets at
                for tgno = 1:length(neuron(tpnno).targets) % for each target
                    switch neuron(tpnno).targets(tgno).to_ntype{:}
                        case 'II'
                            newspike = [ts + neuron(tpnno).targets(tgno).delaysamps neuron(tpnno).targets(tgno).to_synno] ;
                            % insert newspike into IIneuron(neuron(tpnno).targets(tgno).to_nno).inputs
                            IIneuron(neuron(tpnno).targets(tgno).to_nno).inputs = [newspike; IIneuron(neuron(tpnno).targets(tgno).to_nno).inputs] ;
                            IIneuron(neuron(tpnno).targets(tgno).to_nno).inputs = sortrows(IIneuron(neuron(tpnno).targets(tgno).to_nno).inputs) ;
                        case 'TPN'
                            newspike = [ts + neuron(tpnno).targets(tgno).delaysamps neuron(tpnno).targets(tgno).to_synno] ;
                            % apical or basal? Additive or shunting?
                            switch neuron(tpnno).targets(tgno).to_syntype{:}
                                case 'A'
                                    apical(neuron(tpnno).targets(tgno).to_nno).apicalinputs = [newspike; apical(neuron(tpnno).targets(tgno).to_nno).apicalinputs] ;
                                    apical(neuron(tpnno).targets(tgno).to_nno).apicalinputs = sortrows(apical(neuron(tpnno).targets(tgno).to_nno).apicalinputs,1) ;
                                case 'B'
                                    basal(neuron(tpnno).targets(tgno).to_nno).basalinputs = [newspike; basal(neuron(tpnno).targets(tgno).to_nno).basalinputs] ;
                                    basal(neuron(tpnno).targets(tgno).to_nno).basalinputs = sortrows(basal(neuron(tpnno).targets(tgno).to_nno).basalinputs,1) ;
                                case 'AS'
                                    shunts(neuron(tpnno).targets(tgno).to_nno).apicalshuntinputs = [newspike; shunts(neuron(tpnno).targets(tgno).to_nno).apicalshuntinputs] ;
                                    shunts(neuron(tpnno).targets(tgno).to_nno).apicalshuntinputs = sortrows(shunts(neuron(tpnno).targets(tgno).to_nno).apicalshuntinputs, 1) ;
                                case 'BS'
                                    shunts(neuron(tpnno).targets(tgno).to_nno).basalshuntinputs = [newspike; shunts(neuron(tpnno).targets(tgno).to_nno).basalshuntinputs] ;
                                    shunts(neuron(tpnno).targets(tgno).to_nno).basalshuntinputs = sortrows(shunts(neuron(tpnno).targets(tgno).to_nno).basalshuntinputs, 1) ;
                            end % switch
                        otherwise
                    end % switch
                end % for
            end % isfield
        end
    end
    % % and do the same for the II neurons.
    for IIno = 1:simulation.N_IIs
        [IIspike, IIneuron] = II_runstep(ts, IIno, IIneuron, simulation) ;
        if (isfield(IIneuron(IIno), "targets"))
            if IIspike
                % process spike by supplying spikes to IIneuron(IIno).targets at
                % delay IIneuron(IIno).targets().delay
                % process spike by supplying spikes to neuron(tpnno).targets at
                for tgno = 1:length(IIneuron(IIno).targets) % for each target
                    switch IIneuron(IIno).targets(tgno).to_ntype{:}
                        case 'II'
                            newspike = [ts + IIneuron(IIno).targets(tgno).delaysamps IIneuron(IIno).targets(tgno).to_synno] ;
                            % insert newspike into IIneuron(IIneuron(IIno)).targets(tgno).to_nno).inputs
                            IIneuron(IIneuron(IIno).targets(tgno).to_nno).inputs = [newspike; IIneuron(IIneuron(IIno).targets(tgno).to_nno).inputs] ;
                            IIneuron(IIneuron(IIno).targets(tgno).to_nno).inputs = sortrows(IIneuron(IIneuron(IIno).targets(tgno).to_nno).inputs) ;
                        case 'TPN'
                            newspike = [ts + IIneuron(IIno).targets(tgno).delaysamps IIneuron(IIno).targets(tgno).to_synno] ;
                            % apical or basal? Additive or shunting?
                            switch IIneuron(IIno).targets(tgno).to_syntype{:}
                                case 'A'
                                    apical(IIneuron(IIno).targets(tgno).to_nno).apicalinputs = [newspike; apical(IIneuron(IIno).targets(tgno).to_nno).apicalinputs] ;
                                    apical(IIneuron(IIno).targets(tgno).to_nno).apicalinputs = sortrows(apical(IIneuron(IIno).targets(tgno).to_nno).apicalinputs,1) ;
                                case 'B'
                                    basal(IIneuron(IIno).targets(tgno).to_nno).basalinputs = [newspike; basal(IIneuron(IIno).targets(tgno).to_nno).basalinputs] ;
                                    basal(IIneuron(IIno).targets(tgno).to_nno).basalinputs = sortrows(basal(IIneuron(IIno).targets(tgno).to_nno).basalinputs,1) ;
                                case 'AS'
                                    shunts(IIneuron(IIno).targets(tgno).to_nno).apicalshuntinputs = [newspike; shunts(IIneuron(IIno).targets(tgno).to_nno).apicalshuntinputs] ;
                                    shunts(IIneuron(IIno).targets(tgno).to_nno).apicalshuntinputs = sortrows(shunts(IIneuron(IIno).targets(tgno).to_nno).apicalshuntinputs, 1) ;
                                case 'BS'
                                    shunts(IIneuron(IIno).targets(tgno).to_nno).basalshuntinputs = [newspike; shunts(IIneuron(IIno).targets(tgno).to_nno).basalshuntinputs] ;
                                    shunts(IIneuron(IIno).targets(tgno).to_nno).basalshuntinputs = sortrows(shunts(IIneuron(IIno).targets(tgno).to_nno).basalshuntinputs, 1) ;
                            end % switch
                        otherwise
                    end % switch
                end % for
            end %if IIspike
        end % if isfield
    end
    % if TPNspike or IIspike then do spike-time based adaptation for this
    % timestep
end

%% simulation complete.
spikelist = createspikelist(simulation, neuron, IIneuron) ;
if spikesave
    save( strcat("spikelist",string(datetime("today")), ".mat") , "spikelist") ;
end
if spikeplot
    % plot spikes
    figure ;
    spikeraster(spikelist)  ;

    set(groot,'defaultLineLineWidth',2.0) ;

    for tpnno = 1:simulation.N_TPNs
        figure ;
        plot(ahactiv(tpnno,:)', 'DisplayName','ahactiv') ;
        hold on
        plot(neuron(tpnno).TP_threshold, 'DisplayName','threshold') ;
        plot(apicalactivation(tpnno,:), 'DisplayName','ap activ') ;
        plot(basalactivation(tpnno,:), 'DisplayName','ba activ') ;
        title(['TPN ', num2str(tpnno), ' neuron axon hillock and threshold']) ;
        legend ;
    end
    for IIno = 1:simulation.N_IIs
        figure;
        plot(IIneuron(IIno).activation, 'DisplayName','activation') ;
        hold on
        plot(IIneuron(IIno).II_threshold, 'DisplayName','threshold') ;
        title(['II ', num2str(IIno),  ' neuron axon hillock and threshold']) ;
        legend ;
    end
end

% write out weights if required
if weightsave
    wstrings = split(weightfile, ".") ;
    writeweightstofile(strcat(wstrings(1),"_out.txt"), simulation, apical, basal, IIneuron) ;
end
% write out network (with delays) if required
if networksave
    nstrings = split(connectionfile, ".") ;
    writenetworktofile(strcat(nstrings(1),"_out.txt"), simulation, inputneuron,  neuron, IIneuron) ;
end
% format long
% max(ahactiv)

end
