function psthinfos = Expt_SqWave_sorted( CType, Qual )
%
% Usage: psths = Expt_SqWave_sorted( CType, <Qual> )
%
% Shows only the GOOD QUALITY PSTHs for whatever cell type called within the data for all dates
% Master List of all Expt Info files 
% parse_sqwave run for all data dates for just specified cell type 

if (nargin < 2) || isempty(Qual)    %Defaults to 
	Qual = 3;                                                                                         
end

% Constants to describe
TBIN = 10.0; % ms -- approximate
NTrep = 60;  % generally 1 Hz repeats

rootdir = '/home/hoperetina/data/';

%TBIN = 10.0; % ms -- approximate
%NTrep = 60;  % generally 1 Hz repeats

%load([rootdir 'ExperimentDemographics.mat']);
InfoFileNames = dir([rootdir 'Experiment_Info_Files']);
Nexpts = 0;
for nn = 1:length(InfoFileNames)
	if length(InfoFileNames(nn).name) > 10
		a = find(InfoFileNames(nn).name == '.');
		if ~isempty(a) && strcmp( InfoFileNames(nn).name(a-7:a-1), '_Sorted' )
			Nexpts = Nexpts + 1;
			ExptList{Nexpts} = InfoFileNames(nn).name;
		end
	end
end
fprintf( '\n%d experiments sorted.\n', Nexpts );
TotCells = 0;

for ExpR = 1:Nexpts

    %filename = [rootdir 'Experiment_Info_Files/ExptInfo150', int2str(Edate(ExpR)), '.mat'];                                           %Loads ExptInfo file
    %info = load(filename);
		info = load([rootdir 'Experiment_Info_Files/' ExptList{ExpR}]);
    
    for i = 1:length(info.CellQual)
        Gcells = find(info.CellQual <= Qual);                                                                                         %Finds all cells of the desired quality
        Tcells = find(info.CellTypesNum == CType);                                                                                    %Finds all cells of the desired type
    end
    
    Loca = intersect(Gcells, Tcells);                                                                                                 %Puts all elements common to Gcells and Tcells in Loca     
		Ncells = length(Loca);
		
		%data.stimtype{TotCells+nn} = stimtype;
		%fprintf( '\nProcessing Expt 150%s: %d cells\n', int2str(info.EDate(ExpR)), Ncells );
		fprintf( '%2d: Processing %s: %d cells\n', ExpR, ExptList{ExpR}, Ncells );
		
		if Ncells > 0
			psthname = sprintf('%s%s/%s', rootdir, info.directory, info.sqwavefilename{1} );                   
			SWraw = load( psthname );
    
			trigs = SWraw.triggerdata.times(2:end);
			% weird first trigger, and then everything fine
			dt = mean(diff(trigs(2:end)));
			bint = round(TBIN/dt)*dt; % choose multiple of frame rate

			Nreps = round(length(trigs)/round(NTrep));
    
      for nn = 1:length(Loca)    
				cc = Loca(nn);
        spks = SWraw.spikedata.spiketimes{cc} - trigs(1);
				
        psthinfos.cellname{TotCells+nn} = sprintf( '%s C%d', info.EDate , Loca(nn) );
				psthinfos.psth{TotCells+nn} = histc(mod(spks,dt*NTrep), 0:bint:NTrep*dt )/Nreps/bint*1000;
            
				%frsth = psths{cc}(1:((length(psths{cc})-1)/2));  %Gives On:Off Ratios for the desired cells in a vector          
				%secdh = psths{cc}(((length(psths{cc})-1)/2):end-1);
				%pON = max(frsth);
				%pOFF = max(secdh);
				%ONness(end+1) = pON/(pON+pOFF);
			end
			TotCells = TotCells + length(Loca);  
		end
end
psthinfos.dt = TBIN/1000;

if TotCells == 0
	data = [];
	disp( 'No cells found meeting the criteria.' )
else
	fprintf( '\n%d cells found of type %d (qual <= %d).\n\n', TotCells, CType, Qual );
end

