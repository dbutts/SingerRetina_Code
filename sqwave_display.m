function sqwave_display( psthinfo, cellstoplot )
%
% Usage: sqwave_display( psthinfo, <cellstoplot> )

if (nargin < 2) || isempty(cellstoplot)
	cellstoplot = 1:length(psthinfo.psth);
end

Nrows = 3;
Ncols = 4;

if isempty(psthinfo)
	disp('No cells.')
	return
end

dt = psthinfo.dt;
%Ncells = length(psthinfo.psth); 
Ncells = length(cellstoplot); 
cellsplotted = 0;
for nn = 1:ceil(Ncells/(Nrows*Ncols))
	figure
	for mm = 1:min([Ncells-cellsplotted (Nrows*Ncols)])  
		cc = cellsplotted+cellstoplot(mm);
		subplot(Nrows,Ncols,mm); hold on
		NT = length(psthinfo.psth{cc});
		plot((0:(NT-1))*dt,psthinfo.psth{cc})
		title( psthinfo.cellname{cc} )
    %axis([0 NTrep*dt 0 max([40 max(psths{cc})])])
		axis tight
		ylim( [0 max([40 max(psthinfo.psth{cc})])] )
		set(gca,'XTick',0:dt:(NT*dt),'XTickLabel',[])
		plot([1 1]*NT*dt/2,[0 max([40 max(psthinfo.psth{cc})])],'r')
	end
	cellsplotted = cellsplotted + (Nrows*Ncols);
end
