function [handle] = new_fig(varargin)
% [handle] = new_fig(position)
%   create a new figure and return the handle
%
% Inputs:
%  (optional) position: [x_pos,y_pos,x_size,y_size] in inches
%  if no position is specfified, will just return the figure handle
%
% Output: figure handle
%
% Eugene Du
% Feb. 18, 2018

 if nargin == 1 
      handle = figure('units','inches','position',varargin{1});
 else
      handle = figure;
 end

end

