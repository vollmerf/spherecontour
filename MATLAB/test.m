% File    : test.m
% System  : Matlab/Octave
% Purpose : test program for spherecontour.m
% Author  : Frederick W. Vollmer
% Date    : 24 Sep 2014
% Update  : 08 Apr 2020
% Notice  : Copyright 2014-2020
% License : See LICENSE
%
% The algorithms used in this software are described in:
%
%  __F.W. Vollmer, 1995. C program for automatic contouring of spherical 
%  orientation data using a modified Kamb method: Computers & Geosciences, 
%  v. 21, n. 1, p. 31-49.__
%  
% which should be cited by publications using this code, algorithm, or 
% derivative works to produce figures or other content. 
%-------------------------------------------------------------------------

% get comma delimited test file
[filename, pathname] = uigetfile( {'*.csv'});
m = csvread([pathname,filename]);
if isempty(strfind(filename, '_sd'))
  opts = 'dec,inc,cint,mud'; % declination, inclination in degrees
else
  opts = 'str,dip,cint,mud'; % strike, dip in degrees
end

% create a Schmidt plot (lower hemisphere equal-area projection)
[points,lines,frame,grid] = spherecontour(m,opts,5,50);

% set up figure
figure;
hold on;
axis([-1.0 1.0 -1.0 1.0]);
axis('equal');
axis('off');

% define colormap 
% problem: Matlab/Octave sets blanked (NaN) values to first colormap element

% solution 1 - choose colormap with first element = white:
%cmap = gray(256);
%colormap(flipud(cmap));
%lc = 256;
%c1 = [1,1,1];
%c2 = [1,0,0];
%clin = [linspace(c1(1),c2(1),lc)', linspace(c1(2),c2(2),lc)', ... 
%        linspace(c1(3),c2(3),lc)'];
%colormap(clin);

% solution 2 - replace NaN with -0.1 and first colormap element=white
%cmap = jet(256);
%cmap(1,:) = [1 1 1];
%colormap(cmap);
%grid(isnan(grid)) = -0.1;

% solution 3 - define WBGYR colormap, image outside plot is white
T = [255,255,255;0,0,255;0,255,0;255,255,0;255,0,0]./255; 
x = [0,64,128,192,255];
cmap = interp1(x/255,T,linspace(0,1,255));
colormap(cmap);

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

% plot frame, returned as array of (x1, y1, x2, y2), first four are ticks
[n,m] = size(frame);
for i = 1:n
  lx = [frame(i,1), frame(i,3)];
  ly = [frame(i,2), frame(i,4)]; 
  if i <= 4 % ticks
    line ('XData', lx, 'YData', ly, 'Color', 'k', 'LineWidth', 2);
  else % circle
    line ('XData', lx, 'YData', ly, 'Color', 'k', 'LineWidth', 2);
  end
end

% plot points, returned as array of (x,y)
px = points(:,1); 
py = points(:,2); 
h = plot(px, py, 'o');
set(h(1),'MarkerEdgeColor','k','MarkerFaceColor','w', 'MarkerSize', 8)

hold off;
