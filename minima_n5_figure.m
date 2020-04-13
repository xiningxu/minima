% Initialization
clear; clc;
x = linspace(0, pi, 801);
y = linspace(0, pi, 801);
[X, Y] = meshgrid(x,y);

% Constraints
f_con = @(s,t)( cos(s) + cos(t) + cos(s-t) ); 
con = ( f_con(X,Y) <= 1/2);
% Value of the Object Function
val = nan*zeros(size(X));
val_temp = fval(X,Y);
val(con) = real( val_temp(con) );

% Figure
fig = figure(1);
axes1 = axes('Parent',fig);
hold on
contourf(X, Y, val, 400, 'Linestyle', 'none');
col_bar = colorbar();
colormap Jet;
axis('square');
box('on')
xlabel('$$\theta_2$$', 'Interpreter','latex','fontsize', 16);
ylabel('$$\theta_5$$', 'Interpreter','latex','fontsize', 16);

tk = (0:0.25:1)*pi;
tkl = {'0', '\pi/4','\pi/2','3\pi/4','\pi'};
set(axes1, 'Color', [0.8,0.8,0.8], 'Linewidth',1.5,...
    'Layer','top', 'TickDir', 'out','fontsize',12,...
    'Xtick',tk, 'Xticklabel', tkl, 'Ytick',tk,'Yticklabel', tkl);
set(col_bar, 'Limits',[14,15.5], 'Linewidth', 1.5);
% Markers
x1 = acos(-1/4); x2 = acos(-7/8);
pt = [x1, 0; 0, x1; x1, x1; x2, 0; 0, x2; x2, x2];
plot(pt(:,1), pt(:,2), 's', 'Markersize', 25, 'Linewidth',2,...
    'MarkerEdgeColor', [1, 0.41, 0.16]);

Pos = [358,290,30,30; 358,65,30,30; 100,320,30,30; 270,205,30,30;...
        270,65,30,30; 100,235,30,30];
for i = 1:length(pt)
    annotation(fig,'textbox','Units','pixels','Position', Pos(i,:),...
        'Color', [1,0.41,0.16], 'FontName','Arial',...
        'FontWeight','bold','FontSize',18,...
        'VerticalAlignment','middle','HorizontalAlignment','center',...
        'String', num2str(i),...
        'FitBoxToText','off', 'EdgeColor','none')
end

%% Object Function
function val = fval(x, y)
f1c = @(s)( 1 - cos(s) );
f2c = @(s,t)( 3 + cos(s) + cos(t) );
f2s = @(s,t)( sin(s) + sin(t) );
f3c = @(s,t)( cos(s) + cos(t) + cos(s-t) );

r2 = 4 ./ ( 3 + 2*f3c(x,y) );
r = sqrt(r2-1);

A1 = sqrt( 2*f1c(x) ) + sqrt( 2*f1c(y) ) + sqrt( 2*f1c(x-y) );
A2 = sqrt( f2c(x, y) - r .* f2s(x, y) )...
    + sqrt( f2c(x, y) + r .* f2s(x, y) )...
    + sqrt( f2c(x, x-y) - r .* f2s(x, x-y) )...
    + sqrt( f2c(x, x-y) + r .* f2s(x, x-y) )...
    + sqrt( f2c(y, y-x) - r .* f2s(y, y-x) )...
    + sqrt( f2c(y, y-x) + r .* f2s(y, y-x) );
A3 = sqrt( 1 - 2*f3c(x,y));
val = A1 + A2 + A3;
end