function DStuning = Expt_DSdata_sorted( CType, Qual, DSbinning )
%
% Usage: DStuning = Expt_DSdata_sorted( CType, <Qual>, <DSbinning (in sec)> )

% Inputs:
%   Ctype : cell type
%   Qual : Quality of cells (1(= best) to 5) 
%   DSbinning: Bin-size to calculate std of firing rate

% DSexpt order: 0 -90 180 90 45 -135 135 -45 (which was mistakenly 135 again before 1/6/2016)
DIRorder =     [1   7   5  3  2    6   4   8];

Ndirs = length(DIRorder);
Nblocks = 300;

if (nargin < 2) || isempty(Qual)                                                                                                      %Defaults to good quality only
  Qual = 3;                                                                                         
end
if nargin < 3
	DSbinning = 1.0; % sec
end

rootdir = '/home/hoperetina/data/';

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
		
		ExptMonth = str2num(info.EDate(1:2));  % assume if after June, then has DS problem
		
		%data.stimtype{TotCells+nn} = stimtype;
		%fprintf( '\nProcessing Expt 150%s: %d cells\n', int2str(info.EDate(ExpR)), Ncells );
		fprintf( '%2d: Processing %s: %d cells\n', ExpR, ExptList{ExpR}, Ncells );
		
		if ~isfield(info,'DSfilenames')
			for nn = 1:Ncells
				DStuning{TotCells+nn}.cellnames = sprintf( '%s C%d', info.EDate , Loca(nn) );
				DStuning{TotCells+nn}.filenames = [];
				DStuning{TotCells+nn}.DStun = [];
				DStuning{TotCells+nn}.DSstd = [];
			end
		else
			Nspeeds = length(info.DSfilenames);
			if Ncells > 0	
				for nn = 1:Ncells
					DStuning{TotCells+nn}.cellnames = sprintf( '%s C%d', info.EDate , Loca(nn) );
					DStuning{TotCells+nn}.filenames{Nspeeds} = [];
					DStuning{TotCells+nn}.DStun = zeros(Ndirs,Nspeeds);
					DStuning{TotCells+nn}.DSstd = zeros(Ndirs,Nspeeds);
				end
				for mm = 1:Nspeeds
					if ~isempty(info.DSfilenames{mm})
						filename = sprintf('%s%s/%s', rootdir, info.directory(1:end), info.DSfilenames{mm});
						data = load( filename );
						DStuning{TotCells+nn}.filenames{mm} = info.DSfilenames{mm}; % could just extract speed...
					
						trgs = data.triggerdata.times;
						assert( length(trgs)==2400, 'Something different about DS-datafile.' ) 
						dt = mean(diff(trgs(1:Nblocks)));
						T = dt*Nblocks/1000;

						for nn = 1:Ncells
							spks = data.spikedata.spiketimes{Loca(nn)};
							for dd = 1:Ndirs
								%Dspks = spks((spks > trgs((dd-1)*Nblocks+1)) & (spks <= (trgs(dd*Nblocks)+dt))) - trgs((dd-1)*Nblocks+1);
								DIRcount = histc( spks, trgs((dd-1)*Nblocks+1):(DSbinning*1000):trgs(dd*Nblocks)+dt );
								if DIRcount(end) == 0
									DIRcount = DIRcount(1:end-1);
								end
								DStuning{TotCells+nn}.DStun(DIRorder(dd),mm) = mean(DIRcount) / DSbinning;
								DStuning{TotCells+nn}.DSstd(DIRorder(dd),mm) = std(DIRcount) / DSbinning;
							end
						end
					end
				end
				% Correct experiments 7/2015-12/2015 for -45 being replaced by second run of 135
				if ExptMonth >= 7
					disp('    Correcting for DS experiment error (-45 -> 135)' )
					for cc = TotCells+(1:Ncells)
						extra135mean = DStuning{cc}.DStun(8,:);  extra135std = DStuning{cc}.DSstd(8,:);
						DStuning{cc}.DStun(8,:) = NaN;           DStuning{cc}.DSstd(8,:) = NaN;
						DStuning{cc}.DStun(4,:) = (DStuning{cc}.DStun(4,:) + extra135mean)/2;
						DStuning{cc}.DSstd(4,:) = (DStuning{cc}.DSstd(4,:) + extra135std)/2;
					end
				end

			end
		end
				
		TotCells = TotCells + Ncells;
end


	
if TotCells == 0
	DStuning = [];
	disp( 'No cells found meeting the criteria.' )
else
	fprintf( '\n%d cells found of type %d (qual <= %d).\n\n', TotCells, CType, Qual );
end


