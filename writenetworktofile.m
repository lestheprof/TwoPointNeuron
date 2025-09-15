function writenetworktofile(fname, simulation, inputneuron, neuron, IIneuron)
%writenetworktofile writes network out to a file
%   inputs are name of file, structures simulation, , inputneuronneuron and IIneuron
% started LSS 15 Sept 2025.
%
% count number of network arcs for initialising arrays

arcsfrominputs = 0;
arcsfromTPNs = 0 ;
arcsfromIIs = 0 ;

for i = 1:simulation.N_Inputs
    arcsfrominputs = arcsfrominputs + length(inputneuron(i).targets) ;
end
for i = 1 :simulation.N_TPNs
   arcsfromTPNs = arcsfromTPNs + length(neuron(i).targets) ;
end
for i = 1 :simulation.N_IIs
    arcsfromIIs = arcsfromIIs + length(IIneuron(i).targets) ;
end
arcs = arcsfrominputs + arcsfromTPNs + arcsfromIIs ;

% from_ntype	from_nno	to_ntype	to_nno	to_syntype	to_synno	delay
% initialise the entities to be in table
from_ntype = strings(arcs, 1) ;
from_nno = zeros(arcs, 1) ;
to_ntype = strings(arcs, 1) ;
to_nno  = zeros(arcs, 1) ;
to_syntype = strings(arcs, 1) ; 
to_synno = zeros(arcs, 1) ;
delay = zeros(arcs, 1) ;
% fill these arrays
arcno = 1 ;
for i = 1:simulation.N_Inputs
    for tno = 1:length(inputneuron(i).targets)
        from_ntype(arcno) = "X";
        from_nno(arcno) = i ;
        to_ntype(arcno) = inputneuron(i).targets(tno).to_ntype ;
        to_nno(arcno) = inputneuron(i).targets(tno).to_nno;
        to_syntype(arcno) = inputneuron(i).targets(tno).to_syntype ;
        to_synno(arcno) = inputneuron(i).targets(tno).to_synno ;
        delay(arcno) =  inputneuron(i).targets(tno).delay ;
        arcno = arcno + 1 ;
    end
end
for i = 1:simulation.N_TPNs
    for tno = 1:length(neuron(i).targets)
        from_ntype(arcno) = "TPN";
        from_nno(arcno) = i ;
        to_ntype(arcno) = neuron(i).targets(tno).to_ntype ;
        to_nno(arcno) = neuron(i).targets(tno).to_nno;
        to_syntype(arcno) = neuron(i).targets(tno).to_syntype ;
        to_synno(arcno) = neuron(i).targets(tno).to_synno ;
        delay(arcno) =  neuron(i).targets(tno).delay ;
        arcno = arcno + 1 ;
    end
end
for i = 1:simulation.N_IIs
    for tno = 1:length(IIneuron(i).targets)
        from_ntype(arcno) = "II";
        from_nno(arcno) = i ;
        to_ntype(arcno) = IIneuron(i).targets(tno).to_ntype ;
        to_nno(arcno) = IIneuron(i).targets(tno).to_nno;
        to_syntype(arcno) = IIneuron(i).targets(tno).to_syntype ;
        to_synno(arcno) = IIneuron(i).targets(tno).to_synno ;
        delay(arcno) =  IIneuron(i).targets(tno).delay ;
        arcno = arcno + 1 ;
    end
end
 % create table and write it
T= table(from_ntype,	from_nno,	to_ntype,	to_nno,	to_syntype,	to_synno,	delay ) ;
writetable(T, fname,"Delimiter", '\t') ;
end