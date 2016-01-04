function valid_ccs = valid_noise_data( data )
%
% Usage: Exptlist = valid_noise_data( Ndata )
%
% Validity check for noise data (stability of firing rates) -- needs to run on LGN

if isempty(data)
	valid_ccs = [];
	return
end

VAL_THRESH = 0.2;  % tolerates 30% deviation in firing rates
stimdir = '~/Data/SingerRetina/HopeMouse/ProcessedStimuli';
% stimdir = '~/Retina/ProcessedStimuli';

dt = data.dt;
Ncells = length(data.cellname);

valid_ccs = [];
for cc = 1:Ncells
	Nblocks = length(data.stiminfo{cc});
	%spks = [];
	fprintf( '%2d: %s, %d blocks: ', cc, data.cellname{cc}, Nblocks )

	FRs = zeros(Nblocks,1);
	for nn = 1:Nblocks
		stimfilename = sprintf('%s/LongStim%d.mat', stimdir, data.stiminfo{cc}{nn} );
		sdata = load(stimfilename);
		%stim = [stim; sdata.stim;];
		NFR = size(sdata.stim,1);

		rspks = data.spks{cc}{nn};
		rspks = rspks((rspks > 0) & (rspks < (NFR*dt)));
		FRs(nn) = length(rspks)/(NFR*dt);
		fprintf( ' %0.2f Hz  ', FRs(nn) );
		%spks = [spks; rspks(:)+T;];
		%T = T + NFR*dt;
	end
	if Nblocks > 1
		devFR = abs(FRs-mean(FRs))/mean(FRs);
		badones = find(devFR > VAL_THRESH);
		if ~isempty(badones) 
			if Nblocks > 2
				% Check to see if first or last is an outliner
				m1 = mean(FRs(1:end-1)); m2 = mean(FRs(2:end));
			end
			for mm = 1:length(badones)
				fprintf( '\n    Block %d bad: dev = %0.2f ', badones(mm), devFR(badones(mm)) )
			end
		else
			valid_ccs(end+1) = cc;
		end
		fprintf('\n')
	end
end

