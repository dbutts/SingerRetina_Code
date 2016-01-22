function ValBlocks = Ndata_valid_blocks( Ndata, adinfo )
%
% Usage: ValBlocks = Ndata_valid_blocks( Ndata, adinfo )
% 
% looks for more than 50% firing rate change between blocks
% -- if three blocks, eliminate outlier (with bias towards first)
% -- if two blocks, eliminate first

THRESH_BAD = 0.5;

if nargin < 2
	adinfo = cell(1,length(Ndata));
else
	if length(adinfo) < length(Ndata)
		adinfo{length(Ndata)} = [];
	end
end

for ctype = 1:length(Ndata)
	if ~isempty(Ndata{ctype})
	fprintf( 'Type = %d:\n', ctype ) 
	for cc = 1:length(Ndata{ctype}.cellname)
		
		Nblocks = length(Ndata{ctype}.stiminfo{cc});
		ValBlocks{ctype}{cc} = 1:Nblocks;
		fprintf( '  Cell %2d: %d blocks.  ', cc, Nblocks );
		
		% Check firing rates
		FRs = zeros(1,Nblocks);
		for nn = 1:Nblocks
			FRs(nn) = length(Ndata{ctype}.spks{cc}{nn})/600;
		end
		
		dFR = abs(diff(FRs)) ./ FRs(1:end-1);
		if sum(dFR > THRESH_BAD) > 0
			fprintf( 'V' )
			if Nblocks == 2
				% then choose second 
				ValBlocks{ctype}{cc} = 2;
			else
				if dFR(1) > THRESH_BAD
					if dFR(2) < THRESH_BAD
						ValBlocks{ctype}{cc} = 2:3;
					else
						% pick max firing rate
						[~,b] = max(FRs(2:3));
						ValBlocks{ctype}{cc} = b+1;
					end
				else
					if dFR(1) < THRESH_BAD
						ValBlocks{ctype}{cc} = 1:2;
					else
						ValBlocks{ctype}{cc} = 2;
					end
				end
			end
			
		end
		% Replace with adinfo if necessary
		if length(adinfo{ctype}) >= cc
			if ~isempty(adinfo{ctype}{cc})
				if adinfo{ctype}{cc} == 0
					ValBlocks{ctype}{cc} = adinfo{ctype}{cc};
					ValBlocks{ctype}{cc} = adinfo{ctype}{cc};
					fprintf( 'R' )
				end
			end
		end	
		disp(sprintf('%d ', ValBlocks{ctype}{cc} ))
	end
	end
end
