function Rdata = Expt_Repdata_sorted( CType, Qual )
%
% Usage: Rdata = Expt_Repdata_sorted( CType, <Qual> )
% Inputs:
% Ctype : cell type
% Qual : Quality of cells (1(= best) to 5) 
if (nargin < 2) || isempty(Qual)                                                                                                      %Defaults to good quality only
    Qual = 3;                                                                                         
end

if (nargin < 2)                                                                                                                       %Defaults to very good quality on-off cells only
    CType = 1;
    Qual = 3;
end

rootdir = '/home/hoperetina/data/';
%rootdir = '/Users/sarvenaz/Documents/MATLAB/codes/HopeRetina_Sorted/';

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
fprintf( '\n%d experiments sorted.', Nexpts );
TotCells = 0;

for ExpR = 1:Nexpts
	
	info = load([rootdir 'Experiment_Info_Files/' ExptList{ExpR}]);
	Gcells = find(info.CellQual <= Qual);   % finds all cells of the desired quality    
	Tcells = find(info.CellTypesNum == CType);  % finds all cells of the desired type
    
	Loca = intersect(Gcells, Tcells);                                                                                                 %Puts all elements common to Gcells and Tcells in Loca     
	Ncells = length(Loca);

	fprintf( '%2d: Processing %s: %d cells\n', ExpR, ExptList{ExpR}, Ncells );
    
	if Ncells > 0		
		for mm = 1:length(info.repfilenames) 
			if ~isempty(info.repfilenames{mm})
				st = 600; % stimulus size for the repeats            
				edata = parse_rep_expt(sprintf('%s%s/%s', rootdir, info.directory(1:end), info.repfilenames{mm}), Loca, st );
				
				for nn = 1:Ncells	
					Rdata.repinfo{TotCells+nn}{mm} = info.reptype(mm);	
					Rdata.spks{TotCells+nn}{mm} = edata.spks{nn};
					Rdata.cellname{TotCells+nn} = sprintf( '%s C%d', info.EDate , Loca(nn) );
				end
			end	
		end
		fprintf('\n');
		Rdata.dt = edata.dt;
	end
	TotCells = TotCells + Ncells;
end


if TotCells == 0
	Rdata = [];
	disp( 'No cells found meeting the criteria.' )
else
	fprintf( '\n%d cells found of type %d (qual <= %d).\n\n', TotCells, CType, Qual );
end


