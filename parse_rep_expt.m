function data = parse_rep_expt( filename, chans, Nframes )
%
% Usage: data = parse_noise_expt( filename, chans, <Nframes> )

if nargin < 3
	Nframes = 600;
end

sNt = load( filename ); % spikes and triggers
trigs = sNt.triggerdata.times;
% trigs1 is off
dt = mean(diff(trigs(2:end)));
NT = length(trigs); % number of triggers

% Check spacing of triggered and outliers
fprintf( '  %d triggers: %0.2f ms average time\n', NT, dt );

data.filename = filename;
data.chans = chans;
data.dt = dt/1000; %sec 

for nn = 1:length(chans)
	spks = (sNt.spikedata.spiketimes{chans(nn)}-trigs(1))/1000; % in sec
	spksN = spks((spks > 0) & (spks < NT*data.dt));
	spksN = mod(spksN, Nframes*(data.dt));
	indx = find((diff(spksN))<0)+1;
	spksN(indx(1:end)) = -1;
	data.spks{1,nn} = [spksN;-1];
end

