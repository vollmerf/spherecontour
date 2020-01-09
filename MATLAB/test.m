% File    : test.m
% System  : Matlab/Octave
% Purpose : test program for spherecontour.m
% Author  : Frederick W. Vollmer
% Date    : 24 Sep 2014
% Update  : 02 Apr 2018 

% get comma delimited test file
[filename, pathname] = uigetfile( {'*.csv'});
m = csvread([pathname, filename]);
if (strcmp(filename, 'Vollmer_1981_m.csv') == 1)
  opts = 'str,dip,deg'; % strike, dip in degrees
else
  opts = 'dec,inc,deg'; % declination, inclination in degrees
end

% create a Schmidt plot (lower hemisphere equal-area projection)
[points, lines, frame, grid] = spherecontour(m, opts, 10, 50);

% set up figure
figure;
hold on;
axis([-1.0 1.0 -1.0 1.0]);
axis('off');

% define colormap 
% problem: Matlab/Octave sets blanked (NaN) values to first colormap element
% solution 1 - choose colormap with first element = white:
%cmap = gray(256);
%colormap(flipud(cmap));
% solution 2 - replace NaN with negative value and set first colormap element = white:
cmap = jet(256);
cmap(1,:) = [1 1 1];
colormap(cmap);
grid(isnan(grid)) = -0.1;

% plot grid as color gradient
imagesc(-1:1, -1:1, grid);

% plot contours, line segments are returned as array of (x1, y1, x2, y2)
[n,m] = size(lines);
for i = 1:n
  lx = [lines(i,1), lines(i,3)];
  ly = [lines(i,2), lines(i,4)]; 
  %h = plot(lx, ly, 'k');
  %set(h, 'LineWidth', 2);
  line ('XData', lx, 'YData', ly, 'Color', 'k', 'LineWidth', 2);
end

% plot frame, returned as array of (x1, y1, x2, y2), first four are tics
[n,m] = size(frame);
for i = 1:n
  lx = [frame(i,1), frame(i,3)];
  ly = [frame(i,2), frame(i,4)]; 
  if i <= 4
    line ('XData', lx, 'YData', ly, 'Color', 'k', 'LineWidth', 2);
  else
    line ('XData', lx, 'YData', ly, 'Color', 'k', 'LineWidth', 4);
  end
end

% plot points, returned as array of (x,y)
px = points(:,1); 
py = points(:,2); 
h = plot(px, py, 'o');
set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','w')

hold off;