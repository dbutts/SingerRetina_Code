function DSplot( rs, stds, figORsub )
% Usage: DSplot( rs, <stds>, <figure OR subplot> )
%  or    DSplot( DSstruct, <figure OR subplot> )
%
% Make many columns to do multiple

clrs = 'bgcrmk';

if length(rs) == 1
	if nargin < 2
		figORsub = 1;
	else
		figORsub = stds;
	end
	DSstruct = rs;
	rs = DSstruct.DStun;
	stds = DSstruct.DSstd;
else
	if nargin < 2
		stds = [];
	end
	if nargin < 3
		figORsub = 1;
	end
end

if length(figORsub) > 1
	subplot(figORsub(1),figORsub(2),figORsub(3))
else
	figure	
end

[ND,Nplots] = size(rs);
if ND == 1
	rs = rs(:);  stds = stds(:);
	ND = length(rs);
end

rs(end+1,:) = rs(1,:);
if ~isempty(stds)
	stds(end+1,:) = stds(1,:);
end
angs = (0:ND)'/ND*360;

if isempty(stds)
	maxamp = max(rs(:));
else
	maxamp = max(rs(:)+stds(:));
end
%plot(maxamp*[-1 1],[0 0],'k')
%plot([0 0],maxamp*[-1 1],'k')
%polar(angs,ones(ND+1,1)*maxamp,'w')
hold on

for nn = 1:Nplots
	plot(angs,rs(:,nn),clrs(nn),'LineWidth',1.5)
	if ~isempty(stds)
		minstd = rs(:,nn)-stds(:,nn);
		minstd(minstd < 0) = 0;
		plot(angs,minstd,sprintf('%c--',clrs(nn)),'LineWidth',0.5)
		plot(angs,rs(:,nn)+stds(:,nn),sprintf('%c--',clrs(nn)),'LineWidth',0.5)
	end
end
set(gca,'XTick',0:45:360,'XTickLabel',{'0',[],'90',[]','180',[],'270',[],'360'})
axis([0 360 0 maxamp*1.05])
