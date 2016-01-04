function [stim,REPspks,dt] = format_repeat_data( data, cc )
%
% Usage: [Rstim,REPspks,dt] = format_repeat_data( data, cc )

stimdir = '~/Data/SingerRetina/HopeMouse/ProcessedStimuli';

dt = data.dt;
Nstims = length(data.repinfo{cc});
if Nstims == 0
	fprintf( 'No repeat data for cell %d\n', cc );
	stim = [];
	spks = [];
	return;
end

fprintf( '%s: %d blocks.\n', data.cellname{cc}, Nstims );
spks = [];
stim = [];

REPspks = [];

for nn = 1:Nstims
	repfilename = sprintf('%s/RepStim%d.mat', stimdir, data.repinfo{cc}{nn} );
	sdata = load(repfilename);
	stim = [stim; sdata.Rstim;];
	NFR = size(sdata.Rstim,1);
	Trep = (NFR)*dt;  

	%% Calculate number of repeats
	spks = data.spks{cc}{nn};
	%NREPS = floor(ceil(max(spks)/Trep)/2)*2;  % assumes factor of 2 in rep number
	NREPS = length(find(spks < 0));
	fprintf( '  %s: %d reps.\n', repfilename, NREPS );
%	for rep = 1:NREPS
%		T0 = (rep-1)*Trep;
%		rspks = spks((spks > T0) & (spks < (T0+Trep))) - T0;
%		REPspks = [REPspks(:); rspks(:); -1;];
%	end
	REPspks = [REPspks(:); spks(:);];
end
