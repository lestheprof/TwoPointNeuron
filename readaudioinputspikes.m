function [traindata, testdata] = readaudioinputspikes(drivinginputinfo)
%readaudioinputspikes reads two sets of AN spikes, one for training and one
%for test
% reads the driving inputs (test and train) in to a structure, and puts
% the appropriate output into testtraininputs.train/testdata.value.
% NB: the output is the first character of the file name.
[~, testraininfo] = readnetwork(drivinginputinfo) ;
ntrain = 1 ;
ntest = 1 ;
% read input_filelist to get the list of files to be processed
% training files first
inputfid = fopen(join(string([testraininfo(1,1).BDtrain testraininfo(1,2).trainfiles]), '/')) ;
fline = fgetl(inputfid) ;
while ischar(fline)
    trainfilelist{ntrain} = fline ;
    fline = fgetl(inputfid) ;
    ntrain = ntrain + 1 ;
end
ntrain = ntrain - 1 ;

% now load the file contents into structure trainindata 

for fileno = 1:ntrain
    % do we want the whole of the structure?
    t1 = load(join(string([testraininfo(1,1).BDtrain trainfilelist{fileno}]), '/'), 'AN' ) ;
    traindata(fileno).AN = t1.AN.ANParams ;
    traindata(fileno).spikesignal = t1.AN.signal ;
    % and add the targets as well
    traindata(fileno).outputvalue = str2num(t1.AN.ANParams.experiment(1)) ;
end

% testing files
inputfid = fopen(join(string([testraininfo(1,3).BDtest testraininfo(1,4).testfiles]), '/')) ;
fline = fgetl(inputfid) ;
while ischar(fline)
    testfilelist{ntest} = fline ;
    fline = fgetl(inputfid) ;
    ntest = ntest + 1 ;
end
ntest = ntest - 1 ;

% now load the testfile structure
for fileno = 1:ntest
    t2 = load(join(string([testraininfo(1,3).BDtest testfilelist{fileno}]), '/'), 'AN' ) ;
    testdata(fileno).AN = t2.AN.ANParams ;
    testdata(fileno).spikesignal = t2.AN.signal ;
    testdata(fileno).outputvalue = str2num(t2.AN.ANParams.experiment(1)) ;
end 


end
