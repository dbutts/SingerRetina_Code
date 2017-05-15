function RGCinformation_display( CLDdata, fits, ee, cc, force_best )
% 
% Usage: RGCinformation_display( CLDdata, fits, ee, cc, force_best )

maxlag = 30;
LNfit = fits{ee}{cc}.LN;
% Determine best model (LN (sup) or ON-OFF)
[a,b] = max(fits{ee}{cc}.LLs(1:3));
if nargin > 4
	b = force_best;
end
switch b
	case 1, Bfit = LNfit; 
	case 2, Bfit = fits{ee}{cc}.OO; 
	case 3, Bfit = fits{ee}{cc}.sup1; 
end
[Rstim,Rspks,DTstim] = format_repeat_data( CLDdata{ee}, cc );
frac = LNfit.stim_params(1).up_fac;  
tres = DTstim/frac;

figure

subplot(4,5,1:2); hold on
plot( (1:60)/60, CLDdata{ee}.SWdata{cc}(1:60), 'k','LineWidth',1 )
plot( [0.5 0.5], [0 max(CLDdata{ee}.SWdata{cc})*1.05],'r', 'LineWidth',0.5 )
axis tight
title( sprintf( 'e%dc%d class %d (%d): LL = %0.3f, %0.3f', ee, cc, fits{ee}{cc}.class(1), fits{ee}{cc}.class(2), fits{ee}{cc}.LLs(1), a ) )

subplot(4,5,3:4); hold on
nLags = LNfit.stim_params(1).dims(1);
plot([-nLags+1 0],[0 0],'k')
plot( 0:-1:(1-nLags), LNfit.subunits(1).kt,'r','LineWidth',1)
axis([-maxlag 0 1.1*min(LNfit.subunits(1).kt) 1.1*max(LNfit.subunits(1).kt)])
set(gca,'YTick',[])
title('LN')

% Calculate repeat predictions
Nstim = 600;
Rmodstim = Rstim(1:Nstim,fits{ee}{cc}.SPindx);
Nreps = length(find(Rspks == -1));

% Observed response
psth = histc( Rspks, (0:frac*Nstim)*tres)/Nreps;
% Predictions from LN and On-Off model
[~,LNpred,LNmi] = LNfit.eval_model( [], Rmodstim );
[~,OOpred,OOmi] = Bfit.eval_model( [], Rmodstim );
PPs = zeros(1,3);
PPs(1) = PredPower( Rspks, LNpred, tres );
PPs(2) = PredPower( Rspks, OOpred, tres );
[a,b] = max(fits{ee}{cc}.LLs);
BBfit = Bfit;
if b == 4
	BBfit = fits{ee}{cc}.sup3;
elseif b == 5
	BBfit = fits{ee}{cc}.sup4;
end
[~,Bpred] = BBfit.eval_model( [], Rmodstim );
PPs(3) = PredPower( Rspks, Bpred, tres );
subplot(4,5,5)
bar(PPs)
title(sprintf('last = %d', b ));

% Display model parameters
subplot(4,5,6); hold on
plot([-nLags+1 0],[0 0],'k')
plot( 0:-1:(1-nLags), Bfit.subunits(1).kt,'b','LineWidth',1)
if length(Bfit.subunits) > 1
	plot( 0:-1:(1-nLags), Bfit.subunits(2).kt,'g','LineWidth',1)
end
axis tight
xlim([-maxlag 0])
set(gca,'YTick',[])

subplot(4,5,7); colormap gray
imagesc(reshape(Bfit.subunits(1).ksp,21,21),[-1 1]*max(abs(Bfit.subunits(1).ksp)))
set(gca,'XTick',[],'YTick',[])
axis square
title('Exc')

% Nonlinearity
subplot(4,5,9); hold on
meangs = mean(OOmi.gint);
stdgs = std(OOmi.gint);
gs = meangs(1)+3*stdgs(1)*(-1:0.1:1);
maxY = gs(end); %(gs(end)-Bfit.subunits(1).NLoffset);
histgs = histc(OOmi.gint(:,1),gs);
plot( gs(1:end-1)+diff(gs(1:2))/2, histgs(1:end-1)/max(histgs)*0.8*maxY, 'c' );
%plot([gs(1) Bfit.subunits(1).NLoffset gs(end)],[0 0 maxY],'b','LineWidth',1)
if gs(1) < 0
	plot([gs(1) 0 gs(end)],[0 0 maxY],'b','LineWidth',1) % filter offsets included in gints, relative to zero
else
	plot([gs(1) gs(end)],[gs(1) maxY],'b','LineWidth',1) % filter offsets included in gints, relative to zero
end
axis tight

if length(Bfit.subunits) > 1
	subplot(4,5,8); colormap gray
	imagesc(reshape(Bfit.subunits(2).ksp,21,21),[-1 1]*max(abs(Bfit.subunits(2).ksp)))
	set(gca,'XTick',[],'YTick',[])
	axis square
	if Bfit.subunits(2).weight > 0
		title('Exc')
	else
		title('Sup')
	end
	
	subplot(4,5,10); hold on
	gs = meangs(2)+3*stdgs(2)*(-1:0.1:1);
	maxY = gs(end); %(gs(end)-Bfit.subunits(2).NLoffset);
	histgs = histc(OOmi.gint(:,2),gs);
	plot( gs(1:end-1)+diff(gs(1:2))/2, histgs(1:end-1)/max(histgs)*0.8*maxY, 'c' );
	%plot([gs(1) Bfit.subunits(2).NLoffset gs(end)],[0 0 maxY],'g','LineWidth',1)
	if gs(1) < 0
		plot([gs(1) 0 gs(end)],[0 0 maxY],'g','LineWidth',1)
	else
		plot([gs(1) gs(end)],[gs(1) maxY],'g','LineWidth',1)
	end
	axis tight	
end

Nstim = 600;
trange = [1 5];
FRrange = max([floor(trange(1)/tres) 1]):min([ceil(trange(2)/tres) Nstim*frac]);

subplot(4,5,11:15); hold on
plot(FRrange*tres,psth(FRrange)/tres,'k','LineWidth',1)
plot(FRrange*tres,LNpred(FRrange)/tres,'r','LineWidth',1)
plot(FRrange*tres,OOpred(FRrange)/tres,'b','LineWidth',1)
axis([trange(1) trange(2) 0 max(psth)/tres*1.05])

if length(Bfit.subunits) > 1
	subplot(4,5,16:20); hold on
	plot([trange(1) trange(2)],[0 0], 'k' )
	plot(FRrange*tres,LNmi.fgint(FRrange,1),'r')
	plot(FRrange*tres,OOmi.fgint(FRrange,1),'b')
	plot(FRrange*tres,OOmi.fgint(FRrange,2),'g')
	axis tight
	xlim([trange(1) trange(2)])
	set(gca,'YTick',[])
end


