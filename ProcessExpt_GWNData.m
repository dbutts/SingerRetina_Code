function data = ProcessExpt_GWNData( ExptInfo )
%
% Usage: data = ProcessExpt_GWNData( ExptInfo )
%
% Detailed explanation goes here


rootdir = '/home/hoperetina/data/';
% rootdir = '/home/dbutts/Retina/'

info = ExptInfo;
Ncells = length(info.CellQual);
assert( Ncells > 0, 'No cells in experiment.' )
assert( isfield(info,'GWNstimtype'), 'No GWN listed in Expt-Info file.' );

fprintf( 'Expt: %s:\n  Processing %d cells, blocks = ', info.directory, Ncells );

if isfield( info, 'Nblocks_to_use' )
	blocks = info.Nblocks_to_use;
	fprintf( '*' );  % shows there was block info
else
	blocks = 1:length(info.noisefilename);
end
disp(sprintf('%d ', blocks ))

data.expt = info.directory;

for mm = blocks
	fprintf('\nBlock %d:\n', mm )
	if ~isempty(info.GWNfilename{mm})
		edata = parse_noise_expt(sprintf('%s%s/%s', rootdir, info.directory(1:end), info.GWNfilename{mm}), 1:Ncells );
		for nn = 1:Ncells
			data.cellname{nn} = sprintf( '%s C%d', info.EDate , nn );
			data.GWNstiminfo{nn}{mm} = info.GWNstimtype(mm);
			data.spks{nn}{mm} = edata.spks{nn};  
		end
	end
end
fprintf('\n');
data.cell_types = info.CellTypesNum;
data.sort_qual = info.CellQual;
data.blocks = blocks;
data.dt = edata.dt;	
