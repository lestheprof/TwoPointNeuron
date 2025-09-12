function  writeweightstofile(fname, simulation, apical, basal, IIneuron)
%writeweightstofile write(modified) weigts for fine fname so they can be
%re-used

% started LSS 12/09/2025

% Top line: neuron_type	neuron_number	syn_type	syn_number	weight
% and these have to be the names of the arrays to be put into the table
%
% how many weights are there?
nweights = 0 ;
for ntp = 1:simulation.N_TPNs
    nweights = nweights + length(apical(ntp).apicalsynapseweights) ;
    nweights = nweights + length(basal(ntp).basalsynapseweights) ;
end
TPNweights = nweights ;
for iino = 1:simulation.N_IIs
    nweights = nweights + length(IIneuron(iino).weights) ;
end
% initialise the entities to be in table
neuron_type = strings(nweights, 1) ;
neuron_number = zeros(nweights, 1) ;
syn_type = strings(nweights, 1) ; 
syn_number = zeros(nweights, 1) ;
weight = zeros(nweights, 1) ;

% fill arrays
for i = 1:TPNweights
    neuron_type(i) = "TPN";
end
for i= TPNweights + 1 : nweights
    neuron_type(i) = "II" ;
end
wtno = 1 ;
for ntp = 1:length(apical)
    for synno = 1:length(apical(ntp).apicalsynapseweights)
        neuron_number(wtno) = ntp ;
        syn_type(wtno) = "A" ;
        syn_number(wtno) = synno ;
        weight(wtno) = apical(ntp).apicalsynapseweights(synno) ;
        wtno = wtno + 1 ;
    end
end
for ntp = 1:length(basal)
    for synno = 1:length(basal(ntp).basalsynapseweights)
        neuron_number(wtno) = ntp ;
        syn_type(wtno) = "B" ;
        syn_number(wtno) = synno ;
        weight(wtno) = basal(ntp).basalsynapseweights(synno) ;
        wtno = wtno + 1 ;
    end
end
for iino = 1:length(IIneuron)
    for synno = 1:length(IIneuron(iino).weights)
        neuron_number(wtno) = iino ;
        syn_type(wtno) = "S" ;
        syn_number(wtno) = synno ;
        weight(wtno) = IIneuron(iino).weights(synno) ;
        wtno = wtno + 1 ;
    end
end
 % create table and write it
T= table(neuron_type, neuron_number, syn_type, syn_number, weight ) ;
writetable(T, fname,"Delimiter", '\t') ;
end