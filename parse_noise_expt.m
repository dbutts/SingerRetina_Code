function data = parse_noise_expt( filename, chans )
%
% Usage: data = parse_noise_expt( filename, <chans> )

if nargin < 2 
	chans = 1:60;
end

sNt = load( filename );

trigs = sNt.triggerdata.times;
% trigs1 is off
dt = mean(diff(trigs(2:end)));

NT = length(trigs);

% Check spacing of triggered and outliers
fprintf( '\n%d triggers: %0.2f ms average time\n', NT, dt );

data.filename = filename;
data.chans = chans;
data.dt = dt/1000;
for nn = 1:length(chans)
	spks = (sNt.spikedata.spiketimes{chans(nn)}-trigs(1))/1000;
	data.spks{nn} = spks((spks > 0) & (spks < NT*data.dt));
	fprintf( '  %2d (ch %2d): %0.2f Hz\n', nn, chans(nn), length(data.spks{nn})/(NT*dt/1000) );
end

