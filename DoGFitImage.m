function [SpPar, Ksm, Gparams] = DoGFitImage( Ksp, NyNx, Ncov, plotresults, initialseg )
%
% Usage: [SpPar, Ksm] = DoGFitImage( Ksp, NyNx, <Ncov>, <plotresults>, <initialseg> )
% 
% Fits circular Gaussian (or difference-of-Gaussians) using mean-squared error to the image:
%  g = exp( -(P1*(x-x0)^2) + P2*(y-y0)^2) + P3*xy )
%  ** Note standard dev of gaussian can be calculated as 1/sqrt(P3/2)
%
% Inputs:
%   Ksp: filter to fit as one-D array, which will be reshaped using NyNx
%   NyNx: [NY NX] dimensions of filter
%   Ncov: whether to fit single Gaussian (1) or DoG (2) to the image
%   plotresults: want to plot the results, set to 1
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
Cx = floor((ctr-1)/NY)+1;
Cy = mod(ctr-1,NY)+1;

% Initial conditions
if Ncov == 1 
	K0 = [Cy Cx A 0.5];
else
	K0 = [Cy Cx sqrt(A) 1 -A/4 2];
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
SpPar.COV = Gparams(4);
SpPar.Amp = Gparams(3);

x = ((0:NX-1)-Gparams(2));
y = ((0:NY-1)-Gparams(1));

if Ncov == 3
	x2 = ((0:NX-1)-Gparams(12));
	y2 = ((0:NY-1)-Gparams(11));
else
	x2 = x;
	y2 = y;
end

Ksm = abs(Gparams(3)) * exp( -( ones(NY,1)*(x.^2)*abs(Gparams(4)) + (y'.^2)*ones(1,NX)*abs(Gparams(4)) ) );

if Ncov > 1
	Ksm = Ksm - abs(Gparams(5)) * exp( -(ones(NY,1)*(x2.^2)/(1/abs(Gparams(4))+abs(Gparams(6))) + (y2'.^2)*ones(1,NX)/(1/abs(Gparams(4))+abs(Gparams(6)))  ));
	SpPar.COV(2) = 1/ (1/abs(Gparams(4))+abs(Gparams(6)));
	SpPar.Amp(2) = abs(Gparams(5));
	if Ncov == 3
		SpPar.cntr(2,:) = Gparams([12 11]);
	else
		SpPar.cntr(2,:) = SpPar.cntr(1,:);	
	end
end

SpPar.stdevs = 1./sqrt(2*SpPar.COV);

Ksm = Ksm(:);

if plotresults

	
	FRAC_HEIGHT = 1/3; % ellipse at what level of Gaussian?
	thetas = (0:64)/(64)*2*pi;

	xs = SpPar.cntr(1)+1 + sqrt(-2*log(FRAC_HEIGHT))*SpPar.stdevs*sin(thetas);
	ys = SpPar.cntr(2)+1 + sqrt(-2*log(FRAC_HEIGHT))*SpPar.stdevs*cos(thetas);
	
% 	
% 	figure; 
% 	subplot(2,Ncol,1); colormap gray
% 	imagesc(reshape(Ksp,NY,NX)/max(abs(Ksp))',[-1 1])
% 	title('Initial')
% 	hold on
% 	plot(xs,ys,'r','LineWidth',0.5);
% 	plot(SpPar.cntr(1)+1,SpPar.cntr(2)+1,'rx')
% 	if NY == NX, axis square; end

	
	
	figure; 
	subplot(2,2,1); colormap gray
	imagesc(reshape(Ksp,NY,NX)/max(abs(Ksp))',[-1 1])
	hold on
	plot(xs,ys,'r','LineWidth',0.5);
	plot(SpPar.cntr(1)+1,SpPar.cntr(2)+1,'rx')
	if NY == NX, axis square; end
	axis tight
	title('Initial')
	
	subplot(2,2,2); colormap gray
	imagesc(reshape(Ksm,NY,NX)/max(abs(Ksm))',[-1 1])
	hold on
	plot(xs,ys,'r','LineWidth',0.5);
	plot(SpPar.cntr(1)+1,SpPar.cntr(2)+1,'rx')
	if NY == NX, axis square; end
	axis tight
	title('Smoothed')

	subplot(2,2,3)
	plot(reshape(Ksp,NY,NX)')
	xlim([1 NX])
	subplot(2,2,4)
	plot(reshape(Ksm,NY,NX)')
	xlim([1 NX])
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
	
	gaussian = abs(K(3)) * exp( -( ones(NY,1)*(x.^2)*abs(K(4)) + (y'.^2)*ones(1,NX)*abs(K(4)) ) );

	if length(K) > 6
		gaussian = gaussian - abs(K(5)) * exp( -(ones(NY,1)*(x2.^2)/(1/abs(K(4))+abs(K(6))) + (y2'.^2)*ones(1,NX)/(1/abs(K(4))+abs(K(6)))) );
	end
	
	MSE = mean((gaussian(:)-Ksp0).^2);

end
