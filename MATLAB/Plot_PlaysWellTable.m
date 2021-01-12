if ~exist('Delta_fun', 'file')
    addpath([pwd, filesep, 'SK_utils'])
end

global p PlaysWellTruthTable

p = 3;

SK_QAOA_p_objfun_helper;

%% plot truth table
set(0,'DefaultAxesFontSize', 14);
figure(1);
clf

[XX, YY] = meshgrid((0:2^p)+0.5, (0:2^p)+0.5);

pwT = rot90(PlaysWellTruthTable,2);

pcolor(XX, YY,blkdiag(pwT,0))

colormap([1,1,1;0,0,0])
axis square
axis ij

if p <= 4
    set(gca,'xtick',1:2^p,'ytick',1:2^p)
else
    temp = 0:2^(p-4):2^p;
    set(gca,'xtick',temp,'ytick',temp)
end


ax1 = gca; % current axes
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none','YTick',[],'TickLength',[.00014 .0014]);

xx = [2.^(p-(1:p)),1];
xticklocations = xx(1)/2+0.5 + cumsum([0,xx(1:end-1)/2+xx(2:end)/2]);
set(ax2,'xtick',xticklocations,'xticklabels',arrayfun(@(p)sprintf('A_%d',p), 1:p+1,'UniformOutput',0),'xlim',[0.5,2^p+0.5])
set(ax2,'ytick',[],'ydir','reverse','ylim',[0.5,2^p+0.5])

hold on
plot([1;1]*cumsum(xx)+0.5,[0,2^p + 0.5],'--r','linewidth',1.5)
plot([0,2^p+0.5],[1;1]*cumsum(xx)+0.5,'--r','linewidth',1.5)
plot([0,2^p+0.5],[0,2^p+0.5],'-m')
hold off

axis square

ax0 = axes('position',[0,0,1,1],'visible','off');
text(0.05,0.94,sprintf('p=%d',p),'fontsize',27,'fontweight','bold');