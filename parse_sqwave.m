function psths = parse_sqwave( filename, NTrep, toplot )
%
% Usage: psths = parse_sqwave( filename, <NTrep>, <toplot> )
%
% toplot is range of cells whose PSTHs to plot. Make empty if not to plot

TBIN = 10.0; % ms -- approximate

if (nargin < 2) || isempty(NTrep)
	NTrep = 60;  % generally 1 Hz repeats
end

if nargin < 3
	toplot = 1:60;
end

SWraw = load( filename );

trigs = SWraw.triggerdata.times(2:end);
% weird first trigger, and then everything fine
dt = mean(diff(trigs(2:end)));
bint = round(TBIN/dt)*dt; % choose multiple of frame rate

Nreps = round(length(trigs)/round(NTrep));

for nn = 1:ceil(length(toplot)/30)
	if ~isempty(toplot)
		figure
	end 
	for mm = 1:min([length(toplot) 30])
		cc = toplot(mm);
		spks = SWraw.spikedata.spiketimes{cc} - trigs(1);
		psths{cc} = histc( mod(spks,dt*NTrep), 0:bint:NTrep*dt )/Nreps/bint*1000;
		if ~isempty(toplot)
			subplot(6,5,mm); hold on
			plot((1:length(psths{cc})-1)*bint/1000,psths{cc}(1:end-1))
			title(sprintf('Cell %d', cc ))
			%axis([0 NTrep*dt 0 max([40 max(psths{cc})])])
			axis tight
			ylim( [0 max([40 max(psths{cc})])] )
			set(gca,'XTick',0:0.1:(NTrep*dt/1000),'XTickLabel',[])
			plot([1 1]*NTrep*dt/2000,[0 max([40 max(psths{cc})])],'r')
		end
	end
	if length(toplot) > 30
		toplot = toplot(31:end);
	end
end

