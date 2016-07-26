function data = ProcessExpt_SqWaveData( ExptInfo, MinSortQual )
%
% Usage: psthinfos = ProcessExpt_SqWaveData( ExptInfo )

if (nargin < 2) || isempty(MinSortQual)  
	MinSortQual = 3;                                                                                         
end

info = ExptInfo;
Gcells = find(info.CellQual <= MinSortQual);   % finds all cells of the desired quality    
Ncells = length(Gcells);
assert( Ncells > 0, 'No cells in experiment met criteria.' )


% Constants to describe
TBIN = 10.0; % ms -- approximate
NTrep = 60;  % generally 1 Hz repeats

rootdir = '/home/hoperetina/data/';

data.expt = info.directory;

fprintf( 'Processing %s: %d cells\n', data.expt, Ncells );

psthname = sprintf('%s%s/%s', rootdir, info.directory, info.sqwavefilename{1} );                   
SWraw = load( psthname );
    
trigs = SWraw.triggerdata.times(2:end);
% weird first trigger, and then everything fine
dt = mean(diff(trigs(2:end)));
bint = round(TBIN/dt)*dt; % choose multiple of frame rate

Nreps = round(length(trigs)/round(NTrep));
    
for nn = 1:length(Gcells)    
	cc = Gcells(nn);
	spks = SWraw.spikedata.spiketimes{cc} - trigs(1);
				
	data.cellname{nn} = sprintf( '%s C%d', info.EDate , Gcells(nn) );
	data.psth{nn} = histc(mod(spks,dt*NTrep), 0:bint:NTrep*dt )/Nreps/bint*1000;
end
data.dt = TBIN/1000;

