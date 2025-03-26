function [IIspike, IIneuron] = II_runstep(ts, IIno, IIneuron, simulation)
% II_runstep run inhibitory interneuron (leaky integrate and fire) interneuron
% for one timestep (ts).
%   Inhibitory interneurons are simple (non-spiny) LIF neurons, with
%   synapses, and axons primarily aimed at TPNs.
%
% Paramaters
% ts timestep number
% IIno number of the II neuron
% IIneuron structure with info on all II neuros
%
% returns IIspike, true if there is a spike at this ts
% IIneuron, , possib;ly updated
% 
% started 28 Nov 2024 LSS
% modified 5 Fen2025 LSS
%
% do we have any incoming spikes at this ts?

while ((IIneuron(IIno).spikeno <= size(IIneuron(IIno).inputs,1)) && ...
        (IIneuron(IIno).inputs(IIneuron(IIno).spikeno,1) == ts)) % calculate for each synapse: is there a new II spike input? 
    % how should we use alphasynapse in II neurons?
    % treat incoming spike * weight * alphafunctioneffect over time as
    % charge transferred in to the single compartment, and charging a
    % capacitor C
    % use alphafunctioneffect to calculate the  voltage, then add this up
    [alphafn, alphalen] = alphafunctioneffect(IIneuron(IIno).alpha_synapse, IIneuron(IIno).C,IIneuron(IIno).weights(IIneuron(IIno).inputs(IIneuron(IIno).spikeno,2)), ...
        IIneuron(IIno).fracleak, simulation.timestep, 3 * length(IIneuron(IIno).alpha_synapse)) ;
    % add the voltage change: too simple, sensitive to simulation.timestep.
    IIneuron(IIno).activation(ts:(ts + alphalen) -1) = ...
        IIneuron(IIno).activation(ts:(ts + alphalen) -1) + (alphafn * ...
        IIneuron(IIno).synapsemultiplier) ;
    IIneuron(IIno).spikeno = IIneuron(IIno).spikeno + 1;
    IIneuron(IIno).alphalen = alphalen ; % to allow use later on in function
end
if (ts > 1)
    % 24 Jan 2025: calculate activation including decrement due to leakage
    decrement = (IIneuron(IIno).activation(ts-1) - IIneuron(IIno).resetvalue) * IIneuron(IIno).fracleak ;
    IIneuron(IIno).activation(ts) = IIneuron(IIno).activation(ts) - decrement ;
end % ts

% decide whether to spike at this timestep, and store spike if so.
if (IIneuron(IIno).activation(ts) > IIneuron(IIno).II_threshold(ts)) % spike!
    IIspike = 1 ;
    IIneuron(IIno).spikecount = IIneuron(IIno).spikecount + 1 ;
    if (IIneuron(IIno).spikecount) > IIneuron(IIno).maxnospikes
        error("II_runstep: maximum spike count in II %d exceeded at timestep %d", IIno, ts) ;
    end
    IIneuron(IIno).spikes(IIneuron(IIno).spikecount) = ts ;
    % update threshold
    IIneuron(IIno).II_threshold(ts:ts + IIneuron(IIno).th_inc_length -1) = IIneuron(IIno).II_threshold(ts:ts + IIneuron(IIno).th_inc_length -1) + ...
        IIneuron(IIno).thresh_increment ;

    % what to do after firing gets inserted here: should the activation get
    % reduced? 
    if (IIneuron(IIno).resetvalue >= -10) % reset value ge -10 imples reset to this value
        IIneuron(IIno).activation(ts:ts+IIneuron(IIno).alphalen) = IIneuron(IIno).resetvalue ;
    else
        % not currently reset at all
    end
    % II neuron spike time adaptation goes here
else
    IIspike = 0 ;
end
end