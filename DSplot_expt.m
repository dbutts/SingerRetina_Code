function DSplot_expt( DSdata, plottype )
%
% Usage: DSplot_expt( DSdata, <plottype> )

if nargin < 2
	plottype = 1;
end

Ncells = length(DSdata);
Nrows = floor(sqrt(Ncells));
Ncols = ceil(Ncells/Nrows);

figure
for nn = 1:Ncells
	if ~isempty(DSdata{nn}.DStun)
		if plottype == 1
			DSplot(DSdata{nn},[Nrows Ncols nn]);
		else
			DSplot2(DSdata{nn},[Nrows Ncols nn]);
		end
		if isfield(DSdata{nn},'type_qual')
			title(sprintf('Type = %d', DSdata{nn}.type_qual(1) ))
		end
	end
end
