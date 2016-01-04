function DStuning = Expt_DSdata_sorted( CType, Qual, DSbinning )
%
% Usage: DStuning = Expt_DSdata_sorted( CType, <Qual>, <DSbinning (in sec)> )

% Inputs:
%   Ctype : cell type
%   Qual : Quality of cells (1(= best) to 5) 
%   DSbinning: Bin-size to calculate std of firing rate

Ndirs = 8;
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
								DStuning{TotCells+nn}.DStun(dd,mm) = sum(DIRcount) / T;
								DStuning{TotCells+nn}.DSstd(dd,mm) = std(DIRcount) / T;
							end
						end
					end
				end
				%fprintf('\n');
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


