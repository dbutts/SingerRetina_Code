function DStuning = DStuning_expt( info )
%
% Usage: DStuning = DStuning_expt( info )
%
% Calculates DS-tuning for cell cells in particular experiment given experiment-info 'info'

% DSexpt order: 0 -90 180 90 45 -135 135 -45 (which was mistakenly 135 again before 1/6/2016)
DIRorder =     [1   7   5  3  2    6   4   8];

Ndirs = length(DIRorder);
Nblocks = 300;
DSbinning = 1.0; % sec
rootdir = '/home/hoperetina/data/';

Ncells = length(info.CellQual);
fprintf( 'Processing DS on %d cells\n', Ncells );
		
if ~isfield(info,'DSfilenames')
	for nn = 1:Ncells
		DStuning{nn}.cellnames = sprintf( '%s C%d', info.EDate , nn );
		%DStuning{nn}.filenames = info.DSfilenames;
		DStuning{nn}.DStun = [];
		DStuning{nn}.DSstd = [];
	end
else
	Nspeeds = length(info.DSfilenames);
	if Ncells > 0	
		for nn = 1:Ncells
			DStuning{nn}.cellnames = sprintf( '%s C%d', info.EDate , nn );
			DStuning{nn}.type_qual = [info.CellTypesNum(nn) info.CellQual(nn)];
			DStuning{nn}.filenames = info.DSfilenames;
			DStuning{nn}.DStun = zeros(Ndirs,Nspeeds);
			DStuning{nn}.DSstd = zeros(Ndirs,Nspeeds);
		end
		for mm = 1:Nspeeds
			if ~isempty(info.DSfilenames{mm})
				filename = sprintf('%s%s/%s', rootdir, info.directory(1:end), info.DSfilenames{mm});
				data = load( filename );
					
				trgs = data.triggerdata.times;
				assert( length(trgs)==2400, 'Something different about DS-datafile.' ) 
				dt = mean(diff(trgs(1:Nblocks)));
				T = dt*Nblocks/1000;

				for nn = 1:Ncells
					spks = data.spikedata.spiketimes{nn};
					for dd = 1:Ndirs
						DIRcount = histc( spks, trgs((dd-1)*Nblocks+1):(DSbinning*1000):trgs(dd*Nblocks)+dt );
						if DIRcount(end) == 0
							DIRcount = DIRcount(1:end-1);
						end
						%DStuning{nn}.DStun(dd,mm) = sum(DIRcount) / T;
						%DStuning{nn}.DSstd(dd,mm) = std(DIRcount) / T;		
						DStuning{nn}.DStun(DIRorder(dd),mm) = mean(DIRcount) / DSbinning;
						DStuning{nn}.DSstd(DIRorder(dd),mm) = std(DIRcount) / DSbinning;
					end
				end
			end
		end
	end
end

% Correct experiments 7/2015-12/2015 for -45 being replaced by second run of 135
ExptMonth = str2num(info.EDate(1:2));  % assume if after June, then has DS problem
if ExptMonth >= 7
	disp('  Correcting for DS experiment error (-45 -> 135)' )
	for cc = 1:Ncells
		extra135mean = DStuning{cc}.DStun(8,:);  extra135std = DStuning{cc}.DSstd(8,:);
		DStuning{cc}.DStun(8,:) = NaN;           DStuning{cc}.DSstd(8,:) = NaN;
		DStuning{cc}.DStun(4,:) = (DStuning{cc}.DStun(4,:) + extra135mean)/2;
		DStuning{cc}.DSstd(4,:) = (DStuning{cc}.DSstd(4,:) + extra135std)/2;
	end
end
