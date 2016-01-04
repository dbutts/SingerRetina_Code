function DStuning = DStuning_expt( info )
%
% Usage: DStuning = DStuning_expt( info )

Ndirs = 8;
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
						DStuning{nn}.DStun(dd,mm) = sum(DIRcount) / T;
						DStuning{nn}.DSstd(dd,mm) = std(DIRcount) / T;		
					end	
				end
			end
		end
	end
end
