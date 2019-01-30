function [slope] = auto_slope(fig,x,y,style)
% [slope] = auto_slope(fig,x,y,syle)
%   auto generate a slope on a figure
%   currently only available for loglog plots
%   log(y) = log(a) + b*log(x), where b is the slope
%   will use the x and y data to draw a triangle indicating slope
%   below and approximately in the center of the data
%
% Input:
%   fig   : figure handle
%   x     : x data column vector
%   y     : y data column vector
%   style : specifier for number display style eg. '%.2f'
%
% Output:
%   slope : value of slope in loglog scale
%
% Eugene Du
% Jan. 29, 2019

n_data = length(x);

% estimate slope using least squares estimate
a = [ones(n_data,1) log(x)]\log(y);
slope = a(2);

% plot a triangle with the slope
hold(fig.CurrentAxes,'on');

x_range = (x(end) - x(1));
x_s = [x(1) + x_range/4, x(1) + x_range/4 + x_range/2.2 ];
y_s = exp(a(1))*x_s.^(a(2));
y_s = y_s*0.85;
loglog(x_s,y_s,'-k','linewidth',1.5)
loglog(x_s,[y_s(1) y_s(1)],'--k','linewidth',1.5)
loglog([x_s(2) x_s(2)],y_s,'--k','linewidth',1.5)
% above triangle
% loglog([x_s(1) x_s(1)],y_s,'--k','linewidth',1.5)
% loglog(x_s,[y_s(2) y_s(2)],'--k','linewidth',1.5)
text(x(1) + x_range/2.2, 0.85*exp(a(1))*(x(1) + x_range/2.5)^(a(2)),sprintf(style,slope),'fontsize',20)

hold(fig.CurrentAxes,'off');
fig.CurrentAxes; % seems to update the figure
end
