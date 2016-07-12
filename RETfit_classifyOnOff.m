function [newstruct,autoclass] = RETfit_classifyOnOff( fs, Twindow )
%
% Usage: [newstruct,autoclass] = RETfit_classifyOnOff( fitstruct, <Twindow> )
%
% This goes through OOfits (LN fit plus On-Off basic) to sort into ON, OFF, and ON-OFF
% <Twindow> is number of time bins where first maximum must occur in (def = 25)
%
% Output:
%	  autoclass = 2 for On-Off, -1 for Off, and +1 for On

if nargin < 2
	Twindow = 25;
end
newstruct = fs;

% Classify each filter as ON-OFF
Konoff = [0 0];
for nn = 1:2
	kt = fs.OnOff.subunits(nn).kt;
	[MaxAmp,Tmax] = max(kt);
	[MinAmp,Tmin] = min(kt);
	if (Tmax > Tmin)
		if (Tmin < Twindow) && (abs(2*MinAmp) > MaxAmp)
			Konoff(nn) = -1;
		end
	else
		if (Tmax < Twindow) && (MaxAmp > abs(MinAmp/2))
			Konoff(nn) = 1;
		end
	end
end

% Check for ON-OFF
if prod(Konoff) == -1
	% Then on-off
	newstruct.class = 2;
	if max(abs(fs.OnOff.subunits(1).kt)) < max(abs(fs.OnOff.subunits(2).kt)) 
		newstruct.OnOff.subunits = fs.OnOff.subunits([2 1]);
	end
	autoclass = 2;
	return
end

% Otherwise Off or On
newstruct.class = 0;
kt = fs.LN.subunits(1).kt;
[MaxAmp,Tmax] = max(kt);
[MinAmp,Tmin] = min(kt);
if (Tmax > Tmin)
	if (Tmin < Twindow)
		newstruct.class = -1;
	end
else
	if (Tmax < Twindow)
		newstruct.class = 1;
	end
end
autoclass = newstruct.class;
