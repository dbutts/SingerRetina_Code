function data = Expt_NoiseData_sorted( CType, Qual )
%
% Usage: data = Expt_NoiseData_sorted( CType, <Qual> )
%
%   Detailed explanation goes here

if (nargin < 2) || isempty(Qual)                                                                                                      %Defaults to very good quality only
    Qual = 3;                                                                                         
end

if (nargin < 2)                                                                                                                       %Defaults to very good quality on-off cells only
    CType = 1;
    Qual = 3;
end

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
	
	info = load([rootdir 'Experiment_Info_Files/' ExptList{ExpR}]);
  
	Gcells = find(info.CellQual <= Qual); % finds all cells of the desired quality      
	Tcells = find(info.CellTypesNum == CType); % finds all cells of the desired type  
	Loca = intersect(Gcells, Tcells);  % puts all elements common to Gcells and Tcells in Loca     
	Ncells = length(Loca);
	Qs = info.CellQual( Loca );
	
	%data.stimtype{TotCells+nn} = stimtype;
	%fprintf( '\nProcessing Expt 150%s: %d cells\n', int2str(info.EDate(ExpR)), Ncells );
	fprintf( '%2d: Processing %s: %d cells, blocks = ', ExpR, ExptList{ExpR}, Ncells );
	
	if isfield( info, 'Nblocks_to_use' )
		blocks = info.Nblocks_to_use;
		fprintf( '*' );  % shows there was block info
	else
		blocks = 1:length(info.noisefilename);
	end
	disp(sprintf('%d ', blocks ))
	
	if Ncells > 0
		for mm = blocks
			if ~isempty(info.noisefilename{mm})
				%edata = parse_noise_expt([rootdir directory noisefilename(mm)], Loca );
				% extra parsing for now of squarewave:
				edata = parse_noise_expt(sprintf('%s%s/%s', rootdir, info.directory(1:end), info.noisefilename{mm}), Loca );

				for nn = 1:Ncells
					data.cellname{TotCells+nn} = sprintf( '%s C%d', info.EDate , Loca(nn) );
					data.stiminfo{TotCells+nn}{mm} = info.stimtype(mm);
					data.spks{TotCells+nn}{mm} = edata.spks{nn};  % Already specific to Loca
					data.Quals(TotCells+nn) = Qs(nn);
					data.blocks{TotCells+nn} = blocks;
				end
			end
		end	
		fprintf('\n');
		data.dt = edata.dt;
	end
	TotCells = TotCells + Ncells;

end

if TotCells == 0
	data = [];
	disp( 'No cells found meeting the criteria.' )
else
	fprintf( '\n%d cells found of type %d (qual <= %d).\n\n', TotCells, CType, Qual );
end

