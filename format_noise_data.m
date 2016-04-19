function [stim,spks,dt] = format_noise_data( data, cc, blocks )
%
% Usage: [stim,spks,dt] = format_noise_data( data, cc, <blocks_to_include> )
%
% list blocks-to-include

stimdir = '~/Data/SingerRetina/HopeMouse/ProcessedStimuli';
%stimdir = '~/ProcessedStimuli';
%stimdir = '~/Retina/ProcessedStimuli';
%stimdir = '/home/hoperetina/data/ProcessedStimuli';

dt = data.dt;
spks = [];
stim = [];

if nargin < 3
	blocks = [];
end

if isempty(blocks)
	if isfield( data, 'blocks' )
		if iscell(data.blocks)
			blocks = data.blocks{cc};
		else
			blocks = data.blocks;
		end
	else
		blocks = 1:length(data.stiminfo{cc});
	end
end
Nstims = length(blocks);

fprintf( '%s, %d blocks: ', data.cellname{cc}, Nstims )
T = 0;
FRs = [];
for nn = 1:length(blocks)
	if isfield( data, 'GWNstiminfo' )
		% Then GWN stimulus (rather than LongStim)
		stimfilename = sprintf('%s/GWNStim%d.mat', stimdir, data.GWNstiminfo{cc}{nn} );
	else
		% Otherwise standard long stimuli
		stimfilename = sprintf('%s/LongStim%d.mat', stimdir, data.stiminfo{cc}{nn} );
	end
	sdata = load(stimfilename);
	stim = [stim; sdata.stim;];
	NFR = size(sdata.stim,1);

	rspks = data.spks{cc}{nn};
	rspks = rspks((rspks > 0) & (rspks < (NFR*dt)));
	FRs(end+1) = length(rspks)/(size(sdata.stim,1)*dt);
	fprintf( ' %0.2f Hz  ', FRs(end) );
	spks = [spks; rspks(:)+T;];
	T = T + NFR*dt;
end
fprintf('\n')

if (std(FRs) > 0.5*mean(FRs))
	fprintf( '  SIGNIFICANT VARIATION IN FIRING RATES\n' );
end

