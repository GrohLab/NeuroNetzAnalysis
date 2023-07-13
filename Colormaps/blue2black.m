function cm = blue2black(n, varargin)
% Colormap: blues

%-- Parse inputs ---------------------------------------------------------%
if ~exist('n', 'var'); n = []; end
if isempty(n)
   f = get(groot,'CurrentFigure');
   if isempty(f)
      n = size(get(groot,'DefaultFigureColormap'),1);
   else
      n = size(f.Colormap,1);
   end
end
%-------------------------------------------------------------------------%

% Data for colormap:
cm = [
	0       	110      	255
	0           44          102
	0         	0       	0
]/255;

% Modify the colormap by interpolation to match number of waypoints.
cm = tools.interpolate(cm, n, varargin{:});

end