function Rdata = ProcessExpt_RepeatData( ExptInfo, MinSortQual )
%
% Usage: Rdata = ProcessExpt_RepeatData( ExptInfo, <MinSortQual> )

rootdir = '/home/hoperetina/data/';

if nargin < 2
	MinSortQual = 3;
end

info = ExptInfo;
Gcells = find(info.CellQual <= MinSortQual);   % finds all cells of the desired quality    
Ncells = length(Gcells);
assert( Ncells > 0, 'No cells in experiment met criteria.' )

NT = 600; % stimulus size for the repeats            


fprintf( 'Expt: %s:\n  Processing %d cells, blocks = ', info.directory, Ncells );

Rdata.expt = info.directory;

for mm = 1:length(info.repfilenames) 
	fprintf('\nBlock %d:\n', mm )
	if ~isempty(info.repfilenames{mm})
		edata = parse_rep_expt(sprintf('%s%s/%s', rootdir, info.directory(1:end), info.repfilenames{mm}), Gcells, NT );
		for nn = 1:Ncells	
			Rdata.repinfo{nn}{mm} = info.reptype(mm);	
			Rdata.spks{nn}{mm} = edata.spks{nn};
			Rdata.cellname{nn} = sprintf( '%s C%d', info.EDate , Gcells(nn) );
		end
	end
end
Rdata.dt = edata.dt;

