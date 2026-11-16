function fig = plot_scene_geometry(cfg,geom)
%PLOT_SCENE_GEOMETRY 当前仿真的真实几何场景、雷达运动与双点局部放大。
% 独立运行：plot_scene_geometry；无需重新生成回波或执行 BP。
% 雷达图标表示平台参考点；未把静止点散射体画成真实车辆或建筑物。
if nargin<1, cfg=default_config(); end
if nargin<2
    % 最后一次有效发射的时刻；不将末帧尾部空闲时间计入合成孔径。
    timeLast=(cfg.nFrames-1)*cfg.frameTime+ ...
        (cfg.nLoops-1)*cfg.nTx*cfg.slotTime+(cfg.nTx-1)*cfg.slotTime;
else
    timeLast=max(geom.time(:));
end
p0=[0,0]; pEnd=cfg.velocity*timeLast;
points=cfg.targets;
font='Microsoft YaHei'; fonts=listfonts;
if ~any(strcmpi(fonts,font)), font='SimHei'; end
blue=[0.08,0.36,0.70]; orange=[0.91,0.39,0.09]; red=[0.78,0.18,0.24];
grey=[0.43,0.48,0.54]; green=[0.08,0.49,0.39];
fig=figure('Name','仿真场景与雷达运动','Color','w', ...
    'Position',[80,60,1420,950],'Visible',cfg.figureVisible);
layout=tiledlayout(fig,2,2,'TileSpacing','loose','Padding','loose');
title(layout,'TDM-MIMO-SAR 仿真场景：雷达沿轨移动，目标保持静止', ...
    'FontSize',18,'FontWeight','bold');
subtitle(layout,sprintf('77 吉赫兹 | 2 发 4 收 | %d 帧相干积累 | 速度 %.1f 米/秒', ...
    cfg.nFrames,norm(cfg.velocity)),'FontSize',12);

% 左图：全局坐标必须采用等比例，使5米距离与0.4225米孔径关系真实。
ax=nexttile(layout,1,[2,1]); hold(ax,'on');
rectangle(ax,'Position',[cfg.xLimits(1),cfg.yLimits(1),diff(cfg.xLimits),diff(cfg.yLimits)], ...
    'FaceColor',[0.95,0.97,1.00],'EdgeColor',[0.62,0.72,0.85], ...
    'LineStyle',':','LineWidth',1.2);
text(ax,-1.48,8.62,'BP 候选成像范围','FontSize',11,'Color',grey);
for k=1:size(points,1)
    plot(ax,[0,points(k,1)],[0,points(k,2)],'--','Color',[0.76,0.79,0.83],'LineWidth',0.8);
end
plot(ax,[p0(1),pEnd(1)],[p0(2),pEnd(2)],'-','Color',blue,'LineWidth',4);
plot(ax,0,0,'o','MarkerFaceColor',orange,'MarkerEdgeColor','w','MarkerSize',9);
plot(ax,pEnd(1),pEnd(2),'s','MarkerFaceColor',blue,'MarkerEdgeColor','w','MarkerSize',8);
quiver(ax,0.12,0,0,2.25,0,'Color',green,'LineWidth',1.8,'MaxHeadSize',0.23);
text(ax,0.30,2.12,{'侧视参考方向','沿 +y 方向'},'Color',green,'FontSize',11);
quiver(ax,-0.10,-0.35,0.95,0,0,'Color',blue,'LineWidth',1.8,'MaxHeadSize',0.45);
text(ax,1.00,-0.35,'移动方向 +x','Color',blue,'FontSize',11);
text(ax,-0.13,-0.08,{'雷达起点','(0, 0)'},'HorizontalAlignment','right', ...
    'VerticalAlignment','top','FontSize',10,'Color',orange);
text(ax,0.63,0.12,'雷达末位置','FontSize',10,'Color',blue);
scatter(ax,points(:,1),points(:,2),65,red,'filled','MarkerEdgeColor','w');
text(ax,0.18,5.15,{'目标 1、2','见右上放大图'},'Color',red,'FontSize',11);
text(ax,points(3,1)-0.65,points(3,2)+0.18, ...
    sprintf('目标 3\n(%.2f, %.2f) 米',points(3,:)), ...
    'HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',11,'Color',red);
text(ax,points(4,1)+0.75,points(4,2)+0.18, ...
    sprintf('目标 4\n(%.2f, %.2f) 米',points(4,:)), ...
    'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',11,'Color',red);
% 单独标注起点到近邻双点的约5米距离，不将其当作双程路程。
plot(ax,[-1.88,-1.88],[0,5],'-','Color',grey,'LineWidth',1);
plot(ax,[-1.96,-1.80],[0,0],'-','Color',grey);
plot(ax,[-1.96,-1.80],[5,5],'-','Color',grey);
text(ax,-2.02,2.5,'约 5 米','Rotation',90,'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','FontSize',11,'Color',grey);
axis(ax,'equal'); xlim(ax,[-2.35,2.35]); ylim(ax,[-0.75,9.15]); grid(ax,'on'); box(ax,'on');
xlabel(ax,'沿轨位置 x（米）'); ylabel(ax,'侧视位置 y（米）');
title(ax,'全局场景（横纵坐标等比例）','FontSize',14);

% 右上：第一组双点，物理位置来自同一个配置文件。
ax=nexttile(layout,2); hold(ax,'on');
p1=points(1,:); p2=points(2,:); mid=(p1+p2)/2; separation=norm(p2-p1);
for k=1:2
    p=points(k,:); yy=linspace(mid(2)-0.08,mid(2)+0.08,50);
    plot(ax,yy*p(1)/p(2),yy,'--','Color',[0.72,0.76,0.81],'LineWidth',1);
end
plot(ax,p1(1),p1(2),'o','MarkerSize',11,'MarkerFaceColor',orange,'MarkerEdgeColor','w');
plot(ax,p2(1),p2(2),'o','MarkerSize',11,'MarkerFaceColor',red,'MarkerEdgeColor','w');
text(ax,p1(1)-0.008,p1(2)+0.017, ...
    sprintf('目标 1\n(%.4f, %.4f) 米\n初始方位角 -0.5°',p1), ...
    'HorizontalAlignment','right','FontSize',10,'Color',orange);
text(ax,p2(1)+0.008,p2(2)+0.017, ...
    sprintf('目标 2\n(%.4f, %.4f) 米\n初始方位角 +0.5°',p2), ...
    'HorizontalAlignment','left','FontSize',10,'Color',red);
yd=mid(2)-0.033;
plot(ax,[p1(1),p1(1)],[p1(2)-0.005,yd-0.005],':','Color',grey);
plot(ax,[p2(1),p2(1)],[p2(2)-0.005,yd-0.005],':','Color',grey);
quiver(ax,mid(1),yd,p1(1)-mid(1),0,0,'Color',grey,'MaxHeadSize',0.55,'LineWidth',1.2);
quiver(ax,mid(1),yd,p2(1)-mid(1),0,0,'Color',grey,'MaxHeadSize',0.55,'LineWidth',1.2);
text(ax,mid(1),yd-0.009,sprintf('间隔 %.2f 厘米',separation*100), ...
    'HorizontalAlignment','center','VerticalAlignment','top','FontSize',12,'Color',grey);
axis(ax,'equal'); xlim(ax,[mid(1)-0.20,mid(1)+0.20]); ylim(ax,[mid(2)-0.075,mid(2)+0.075]);
grid(ax,'on'); box(ax,'on');
xlabel(ax,'沿轨位置 x（米）'); ylabel(ax,'侧视位置 y（米）');
title(ax,'近邻静止点目标局部放大','FontSize',14);

% 右下：13帧真实发射位置；每帧内较密、帧间留有空闲间隔。
ax=nexttile(layout,4); hold(ax,'on');
for f=1:cfg.nFrames
    tt=(f-1)*cfg.frameTime+(0:cfg.nLoops-1)*cfg.nTx*cfg.slotTime;
    plot(ax,cfg.velocity(1)*tt,cfg.velocity(2)*tt,'.','Color',blue,'MarkerSize',8);
end
plot(ax,p0(1),p0(2),'o','MarkerFaceColor',orange,'MarkerEdgeColor','w','MarkerSize',11);
plot(ax,pEnd(1),pEnd(2),'s','MarkerFaceColor',blue,'MarkerEdgeColor','w','MarkerSize',10);
quiver(ax,0,0.065,pEnd(1),pEnd(2),0,'Color',blue,'LineWidth',1.8,'MaxHeadSize',0.25);
text(ax,pEnd(1)/2,0.082,'雷达沿 +x 方向匀速移动','HorizontalAlignment','center','Color',blue,'FontSize',11);
text(ax,-0.012,-0.023,{'起点','(0, 0)'},'HorizontalAlignment','right','FontSize',10,'Color',orange);
text(ax,pEnd(1)+0.012,-0.023, ...
    sprintf('末位置\n(%.4f, %.1f) 米',pEnd), ...
    'HorizontalAlignment','left','FontSize',10,'Color',blue);
text(ax,pEnd(1)/2,-0.095, ...
    sprintf('有效孔径 %.4f 米；最后发射时刻 %.4f 秒',norm(pEnd),timeLast), ...
    'HorizontalAlignment','center','FontSize',11,'Color',grey);
text(ax,pEnd(1)/2,-0.125,'蓝点为各帧 Tx1 的发射位置，帧间空闲段保留', ...
    'HorizontalAlignment','center','FontSize',10,'Color',grey);
axis(ax,'equal'); xlim(ax,[-0.12,0.65]); ylim(ax,[-0.16,0.14]); grid(ax,'on'); box(ax,'on');
xlabel(ax,'沿轨位置 x（米）'); ylabel(ax,'侧视位置 y（米）');
title(ax,'雷达运动轨迹局部放大','FontSize',14);

set(findall(fig,'-property','FontName'),'FontName',font);
if cfg.saveFigures
    if ~exist(cfg.outputDir,'dir'), mkdir(cfg.outputDir); end
    exportgraphics(fig,fullfile(cfg.outputDir,'07_仿真场景与雷达运动.png'),'Resolution',220);
    savefig(fig,fullfile(cfg.outputDir,'07_仿真场景与雷达运动.fig'));
end
if strcmpi(cfg.figureVisible,'off'), close(fig); end
end
