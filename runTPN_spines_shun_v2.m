% script to run TPN_spines_shun.m
duration = 1;
timestep = 0.0001;
n_apicalinputs = 2;
n_basalinputs =6;
tau_apical = 0.005 ;
tau_basal = 0.006;
C_apical = 10^-7 ;
C_basal = 2 * 10^-7 ;
R_apical = 10^5 ;
R_basal = 10^5 ;
R_synap_dendrite = 5 * 10^5 ;
R_synba_dendrite = 5 * 10^5 ;
C_apical_spine = 10^-7 ;
C_basal_spine = 2 * 10^-7 ;
R_apical_spine = 10^5 ;
R_basal_spine = 10^5 ;
R_synap_spine = 5 * 10^5 ;
R_synba_spine =  5 * 10^5 ;
% basalinputs = [0.51 1; 0.7 2; 0.7 1];
% basalinputs=[] ;
basalinputs = [0.51 1; 0.52 2];
% apicalinputs = [0.48 1];
apicalinputs = [0.49 1; 0.51 2];
apicalsynapseweights = [0.5 0.6];
basalsynapseweights =[0.5 0.9 1.1 0.98 1.1 1.3] ;
% shunting synapse materials
noapicalshunts = 3;
nobasalshunts = 3 ;
apicalshuntweights = [0.9 0.8 0.5];
basalshuntweights = [0.9 0.6 0.8];
apicalshuntduration = 0.01 ;
basalshuntduration = 0.012 ;
apicalshuntinputs = [0.48 2; 0.50 1; 0.6 2; 0.51 3] ;
basalshuntinputs = [0.55 1; 0.54 2; 0.56 3; 0.52 1] ;


% spiking stuff
refractoryperiod = 0.006 ; % 6ms
relrefperiod = 0.04 ;% 40ms
thresh_leap = 8 ;
thresh_decay = 100 ;
thresh_value = 5 ;
maxnospikes = 100 ;

    [aa,  ba, ahactiv, threshold,neuron, spikelist] =  TPN_spines_shun_v2('duration',duration, 'timestep',timestep,  'noapicalinputs',n_apicalinputs, ...
    'nobasalinputs', n_basalinputs, 'apicalinputs', apicalinputs, 'basalinputs',basalinputs,  ...
    'apicalsynapseweights', apicalsynapseweights, 'basalsynapsewights', basalsynapseweights,...
    'tau_apical',tau_apical, 'tau_basal', tau_basal, 'C_apical', C_apical, ...
    'C_basal' , C_basal, 'R_apical', R_apical, 'R_basal', R_basal, 'R_synap_dendrite', ...
    R_synap_dendrite, 'R_synba_dendrite', R_synba_dendrite, ...'
    'C_apical_spine', C_apical_spine, ...
    'C_basal_spine' , C_basal_spine, 'R_apical_spine', R_apical_spine, 'R_basal_spine', R_basal_spine, ...
    'R_synap_spine', R_synap_spine, 'R_synba_spine', R_synba_spine, ...
    'Refractoryperiod', refractoryperiod, 'RelRefPeriod', relrefperiod, 'Thresh_Leap', ...
    thresh_leap, 'Thresh_Decay', thresh_decay, 'Thresh_Value', thresh_value, 'maxnospikes', maxnospikes, ...
    'noapicalshunts', noapicalshunts, 'nobasalshunts' , nobasalshunts, 'apicalshuntweights', apicalshuntweights, ...
    'basalshuntweights', basalshuntweights, 'apicalshuntduration', apicalshuntduration, 'basalshuntduration', basalshuntduration, ...
    'apicalshuntinputs', apicalshuntinputs, 'basalshuntinputs', basalshuntinputs) ;

    % plot spikes
    % figure ; spikeraster(spikelist, 'endtime', 1)  ;
    figure ;
    plot(ahactiv);
    hold on
    plot(threshold) ;
    xlim([0.4/timestep 0.6/timestep]) ;
