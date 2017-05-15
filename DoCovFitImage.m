function [SpPar, Ksm, Gparams] = DoCovFitImage( Ksp, NyNx, Ncov, plotresults, initialseg )
%
% Usage: [SpPar, Ksm, Gparams] = DoCovFitImage( Ksp, NyNx, <Ncov>, <plotresults>, <initialseg> )
% 
% Fits elliptical Gaussian (or difference-of-Gaussians) using mean-squared error to the image:
%  g = exp( -(P1*(x-x0)^2) + P2*(y-y0)^2) + P3*xy )
%
% Inputs:
%   Ksp: filter to fit as one-D array, which will be reshaped using NyNx
%   NyNx: [NY NX] dimensions of filter
%   Ncov: whether to fit single Gaussian (1) or DoG (2) to the image
%   plotresults: want to plot the results, set to 1. Make 2 if want to see fitted comparison.
%   initialseg: initial guess for parameters
% Outputs:
%   SpPar: parameters of fit: cntr [y0 x0], COV [P1 P2 P3], amplitude
%   Ksm: fit images
%   Gparams: parameters above ordered to be input as initial guess

if (nargin < 3) || isempty(Ncov)
	Ncov = 1;
end
if nargin < 4
	plotresults = 0;
end
if nargin < 5
	initialseg = [];
end

NX = NyNx(2);  NY = NyNx(1);

% Minimize mean-squared error with single Gaussian fit
[A,ctr] = max(Ksp);
Cx = floor((ctr-1)/NY); %+1;
Cy = mod(ctr-1,NY); %+1;

% Initial conditions
if Ncov == 1 
	K0 = [Cy Cx A 0.5 0.5 0];
else
	K0 = [Cy Cx sqrt(A) 1 1 0 -A/4 2 2 0];
	if Ncov > 2
		K0(end+(1:2)) = [Cy Cx];
	end
end

if ~isempty(initialseg)
	K0(1:length(initialseg)) = initialseg;
end

%% Set up optimization
outeropt.Display = 'off';
%outeropt.Display = 'iter';
%outeropt.Algorithm = 'quasi-newton';

% Minimization
Gparams = fminsearch( @(K) spatial_fit_internal( K, Ksp, NY, NX ), K0, outeropt );

SpPar.cntr = Gparams([2 1]);
SpPar.COV = Gparams(4:6);
SpPar.Amp = Gparams(3)^2;

x = ((0:NX-1)-Gparams(2));
y = ((0:NY-1)-Gparams(1));

if Ncov == 3
	x2 = ((0:NX-1)-Gparams(12));
	y2 = ((0:NY-1)-Gparams(11));
else
	x2 = x;
	y2 = y;
end

Ksm = abs(Gparams(3)) * exp( -(ones(NY,1)*(x.^2)*abs(Gparams(5)) + (y'.^2)*ones(1,NX)*abs(Gparams(4)) + (y'*x)*Gparams(6)) );

if Ncov > 1
	Ksm = Ksm - abs(Gparams(7)) * exp( -(ones(NY,1)*(x2.^2)/(1/abs(Gparams(5))+abs(Gparams(9))) + (y2'.^2)*ones(1,NX)/(1/abs(Gparams(4))+abs(Gparams(8))) + (y2'*x2)*Gparams(10)) );
	SpPar.COV(2,:) = [1./(1./abs(Gparams(5:6))+abs(Gparams(8:9))) Gparams(10)];
	SpPar.Amp(2,1) = Gparams(7);
	if Ncov == 3
		SpPar.cntr(2,:) = Gparams([12 11]);
	else
		SpPar.cntr(2,:) = SpPar.cntr(1,:);	
	end
end

% Calculate major and minor axis  and real widths
SpPar.stdevs_xy = 1./sqrt(2*SpPar.COV(:,1:2));
C = [SpPar.COV(1,[1 3]); SpPar.COV(1,[3 2]);];
[V,D] = eig(C);
if sum(diag(D) < 0) > 0
	SpPar.stdevs = [0 0];
else
	SpPar.stdevs = 1./sqrt(2*diag(D))';
end
SpPar.major_axis = V(:,1)';
SpPar.angle = 180/pi * atan(V(2,1)/V(1,1));

if Ncov > 1
	C = [SpPar.COV(2,[1 3]); SpPar.COV(2,[3 2]);];
	[V,D] = eig(C);
	SpPar.stdevs(2,:) = 1/sqrt(2*diag(D))';
	SpPar.major_axis(2,:) = V(:,1)';
	SpPar.angle(2) = 180/pi * atan(V(2,1)/V(1,1));
end

Ksm = Ksm(:);

if plotresults
	if plotresults > 1
		Ncol = 2;
	else 
		Ncol = 1;
	end
	
	% Make ellipse: gotta be a better way...
	C = [SpPar.COV(1,[1 3]); SpPar.COV(1,[3 2]);];
	thetas = (0:64)/(64)*2*pi;
	rs = zeros(1,length(thetas));
	for nn = 1:length(thetas)
		%rs(nn) = ( SpPar.COV(1)*(cos(thetas(nn))^2) + 2*SpPar.COV(3)*cos(thetas(nn))*sin(thetas(nn)) + SpPar.COV(2)*(sin(thetas(nn))^2) )^(-0.5);
		v = [cos(thetas(nn)) sin(thetas(nn))];
		rs(nn) = (2*v*C*v')^(-0.5);
	end		
	
	FRAC_HEIGHT = 0.5; % ellipse at what level of Gaussian?
	xs = SpPar.cntr(1)+1 + sqrt(-2*log(FRAC_HEIGHT))*rs.*sin(thetas);
	ys = SpPar.cntr(2)+1 + sqrt(-2*log(FRAC_HEIGHT))*rs.*cos(thetas);
	
	
	figure; 
	subplot(2,Ncol,1); colormap gray
	imagesc(reshape(Ksp,NY,NX)/max(abs(Ksp))',[-1 1])
	title('Initial')
	hold on
	plot(xs,ys,'r','LineWidth',0.5);
	plot(SpPar.cntr(1)+1,SpPar.cntr(2)+1,'rx')
	if NY == NX, axis square; end
	
	subplot(2,Ncol,Ncol+1)
	plot(reshape(Ksp,NY,NX)')
	xlim([1 NX])
	if NY == NX, axis square; end

	if plotresults > 1
		subplot(2,Ncol,2); colormap gray
		imagesc(reshape(Ksm,NY,NX)/max(abs(Ksm))',[-1 1])
		title('Smoothed')
		hold on
		plot(xs,ys,'r','LineWidth',0.5);
		if NY == NX, axis square; end
	
		subplot(2,Ncol,4)
		plot(reshape(Ksm,NY,NX)')
		xlim([1 NX])
		if NY == NX, axis square; end
	end
end

end


function MSE = spatial_fit_internal( K, Ksp0, NY, NX )

	x = ((0:NX-1)-K(2));
	y = ((0:NY-1)-K(1));
	if length(K) > 10
		x2 = ((0:NX-1)-K(12));
		y2 = ((0:NY-1)-K(11));
	else
		x2 = x;
		y2 = y;
	end
	
	gaussian = abs(K(3)) * exp( -(ones(NY,1)*(x.^2)*abs(K(5)) + (y'.^2)*ones(1,NX)*abs(K(4)) + (y'*x)*abs(K(6))) );
	if length(K) > 6
		gaussian = gaussian - abs(K(7)) * exp( -(ones(NY,1)*(x2.^2)/(1/abs(K(5))+abs(K(9))) + (y2'.^2)*ones(1,NX)/(1/abs(K(4))+abs(K(8))) + (y2'*x2)*K(10)) );
	end
	
	MSE = mean((gaussian(:)-Ksp0).^2);

end
