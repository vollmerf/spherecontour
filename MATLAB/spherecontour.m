function [points,lines,frame,grid] = spherecontour( ...
  data,options,nlevels,ngrid,cint,sigma)
% SPHERECONTOUR  Contoured spherical projection of directional data.
%   Spherical coordinates are in a right-handed reference frame X=right 
%   (E), Y=top (N), Z=up. Options are set in the string paramenter 
%   'options'.
%
% Input
%   data    : array of [theta, phi] angles, see 'options' for formats.   
%   options : include in string any non-default options:
%             theta values:
%               ''       = longitude (CCW from X=0 at right)
%               'dec'    = azimuth, trend, declination (CW from Y=N)
%               'str'    = strike (CW from Y=N at top)
%               'dir'    = dip direction (CW from Y=N at top)
%             phi values: 
%               ''       = colatitude or zenith (down from Z)
%               'lat'    = altitude or latitude (up from XY plane)
%               'inc'    = inclination or plunge (down from XY plane)
%               'dip'    = dip of plane (down from XY plane)
%               'nad'    = nadir (up from -Z)
%             angle format: 
%               ''       = degrees
%               'rad'    = radians
%               'grd'    = gradians
%             projection:
%               ''       = equal area (Schmidt plot)
%               'stereo' = stereographic (stereogram)
%             hemisphere:
%               ''       = lower hemisphere 
%               'upper'  = upper hemisphere
%             data type:
%               ''       = axes (undirected)
%               'vec'    = vectors (directed)
%             contouring method:
%               ''       = modified Kamb 
%               'mud'    = modified Kamb multiples of uniform density
%               'sch'    = modified Schmidt (not recommended)
%               'ncon'   = no contouring
%             contour spacing:
%               ''       = equal spaced levels, nlevels=10 gives 9 contour 
%                          lines
%               'cint'   = set contour intervals (cint)
%             contour smoothing:
%               ''       = exponential (recommended)
%               'sma'    = inverse area
%               'sms'    = inverse area squared
%               'nsm'    = none
%             grid interpolation: fixes edge effect due to square grid
%               ''       = 5 parts
%               'gi0'    = off
%               'gi2'    = 2 parts
%               'gi3'    = 3 parts
%               'gi4'    = 4 parts
%               'gi5'    = 5 parts
%               'gi6'    = 6 parts
%               'gi8'    = 8 parts
%               'gi10'   = 10 parts
%             frame:     
%	            ''       = circle and tics             
%               'ntc'    = circle without tics              
%	            'nfr'    = no frame 
%             grid:
%               ''       = grid     
%               'ngd'    = no grid  %             
%   nlevels : number of levels spaced over the distribution, 10 will 
%             give 9 contour lines, default = 10, unused for 'cint' option 
%   ngrid   : number of grid nodes, use 30 for draft plots, 50 or more for 
%             final plots, default = 30
%   cint    : contour interval for 'cint' option, default = 1.0 
%   sigma   : Kamb method sigma in standard deviations, default = 3.0
%
% Output
%   points  : projected data points in unit circle as [x,y]
%   lines   : projected contour line segments in unit circle as 
%             [x1,y1,x2,y2], includes tick marks and a circular frame
%   frame   : tick marks and circle as line segments, first four are ticks
%   grid    : grid for display of color gradient   
% 
% Syntax
%   [points]=spherecontour(m);
%   [points,lines]=spherecontour(m);
%   [points,lines,frame]=spherecontour(m);
%   [points,lines,frame]=spherecontour(m,'str,dip',5,50);
%   [points,lines,frame,grid]=spherecontour(m,'dec,inc,mud,cint',5,50);

  switch nargin
    % data,options,nlevels,ngrid,cint,sigma
    case 1
      options = '';
      nlevels = 10;
      ngrid = 30;
      cint = 1.0;
      sigma = 3.0;
    case 2
      nlevels = 10;
      ngrid = 30;
      cint = 1.0;
      sigma = 3.0;
    case 3
      ngrid = 30;
      cint = 1.0;
      sigma = 3.0;
    case 4
      cint = 1.0;
      sigma = 3.0;
    case 5
      sigma = 3.0;
    case 6
      sigma = sigma;
    otherwise
      return
  end  
  cmin = cint;
  cmax = 1000.0; % not used
  if nargout < 4 % no grid
    opts.grid = 1
  end
  if nargout < 3 % no frame
    opts.frame = 2
  end
  if nargout < 2 % no lines
    opts.method = 2
  end
  if nargout < 1 
    return
  end
  [n,m] = size(data);
  dc = zeros(n,3);
  pts = zeros(n,2);
  opts = getOptions(options);
  if (opts.angfmt == 1) % degrees
    f = pi/180.0;
  elseif (opts.angfmt == 2) % gradians
    f = pi/200.0;
  else % radians
    f = 1.0;
  end
  c = 0;
  for i = 1:n
    t = toSphereTheta(f*data(i,1), opts.thetafmt);
    p = toSpherePhi(f*data(i,2), opts.phifmt);
    d = toDirCos(t, p);
    dc(i,:) = d(:);
    [x, y, visible] = sphereProject(d, opts.proj, opts.hemi, opts.directed);
    if visible
      c = c + 1;
      pts(c,1) = x; 
      pts(c,2) = y;
    end 
  end
  points = pts(1:c,:); % truncate
  if opts.method < 2 % kamb or schmidt
    [grid, lines] = contour(dc, ngrid, sigma, nlevels, cint, cmin, cmax, opts);
    if opts.grid == 0
      grid = processGrid(grid, opts.interp);
    end
    frame = drawFrame(opts.frame);
  end
end

function [grid] = processGrid(grid, interp)
% PROCESSGRID  Process the final grid for output.
% The number of grid nodes is increased and intermediate values 
% interpolated. This is done to fix edge effects at the perimeter due to 
% the square grid. Nodes outside the circle are set to NaN.
  grid = grid';
  [n,m] = size(grid);
  [x y] = meshgrid(1:n);
  if (interp < 0.0)
    zi = grid;
  else
    [xi yi] = meshgrid(1:interp:n);
    zi = interp2(x,y,grid,xi,yi);
  end
  [ni,mi] = size(zi);
  r = 0.5 * (ni-1);
  r2 = r*r;
  [xi yi] = meshgrid(1:ni);
  % NaN nodes outside circle
  zi((xi-r-1.0).^2 + (yi-r-1.0).^2 > r2) = NaN;
  grid = zi;
end

function [x, y, visible] = sphereProject(dc, proj, hemi, directed) 
% SPHEREPROJECT  Projects direction cosines to cartesian coordinates of 
% unit spherical projection.                                             
  t = dc;
  if (hemi == 0) % lower 
    t(3) = -t(3);
  end
  if (directed == false) && (t(3) < 0.0)
    t = -t;
  end
  if (t(3) < 0.0)
    x = 99.0;
    y = 99.0;
    visible = false; 
  else
    if (proj == 1) % stereographic 
      f = 1.0/(1.0+t(3));
    else % equal area
      f = 1.0/sqrt(1.0+t(3));
    end
    x = f*t(1); 
    y = f*t(2);
    visible = true; 
  end
end

function [dc] = sphereBProject(x, y, proj, hemi) 
% SPHEREBPROJECT  Back projects Cartesian coordinates of unit spherical 
% projection to direction cosines.                                       
  r2 = (x*x)+(y*y);
  if (proj == 1) % stereographic  
    dc(3) = (1.0-r2)/(1.0+r2); 
    f = 1.0+dc(3); 
  else % equal area
    f = sqrt(abs(2.0-r2)); 
    dc(3) = 1.0-r2; 
  end
  dc(1) = f*x; 
  dc(2) = f*y;
  if (hemi == 0) % lower 
    dc(3) = -dc(3);
  end
end

function [dc] = toDirCos(theta, phi) 
% TODIRCOS  Converts theta, phi in radians to XYZ direction cosines. 
  s = sin(phi);
  dc(1) = cos(theta) * s;
  dc(2) = sin(theta) * s;
  dc(3) = cos(phi);   
end

function [spherePhi] = toSpherePhi(phi, fmt)
% TOSPHEREPHI  Convert from user coordinates in radians to colatitude 
% (phi) in radians.
  switch fmt 
    case 0 % colatitude
      p = phi;
    case 1 % altitude, latitude
      p = 0.5*pi - phi;
    case 2 % inclination
      p = phi + 0.5*pi;
    case 3 % dip
      p = pi - phi;
    case 4 % nadir
      p = pi - phi;
    otherwise % colatitude
      p = phi;
  end
  spherePhi = p;
end

function [sphereTheta] = toSphereTheta(theta, fmt)
% TOSPHERETHETA  Convert from user coordinates in radians to longitude 
% (theta) in radians.
  switch fmt
    case 0 % longitude
      t = theta;
    case 1 % azimuth
      t = 0.5*pi - theta;
    case 2 % strike
      t = pi - theta;
    case 3 % dip direction
      t = 1.5*pi - theta;
    otherwise 0 % longitude
      t = theta;
  end
  sphereTheta = t;
end

function [t1, t2, visible] = lineCircleInt(x1, y1, x2, y2, xc, yc, r) 
% LNECIRCLEINT  Determine intersection parameters for line segment and 
% circle. Adopted from Rankin 1989, p.220.                               
  visible = false;
  t1 = 0.0;
  t2 = 1.0;
  dx = x2-x1; 
  dy = y2-y1; 
  dxc = x1-xc; 
  dyc = y1-yc;
  a = dx*dxc + dy*dyc; 
  b = dx*dx + dy*dy; 
  c = dxc*dxc + dyc*dyc - r*r;
  disc = a*a - b*c;
  if ((disc > 0.0) && (abs(b) > 1e-9)) 
    d = sqrt(disc);
    t1 = (-a + d)/b; 
    t2 = (-a - d)/b;
    if (t1 > t2) 
      t = t1; 
      t1 = t2; 
      t2 = t; 
    end
    visible = true; 
  end
end

function [cx1, cy1, cx2, cy2, visible] = clipLineCircle(xc, yc, r, x1, y1, x2, y2) 
% CLIPLINECIRCLE  Clip line segment to circle. 
  cx1 = x1;
  cy1 = y1;
  cx2 = x2;
  cy2 = y2;
  visible = false; 
  if (((x1 < xc-r) && (x2 < xc-r)) || ((x1 > xc+r) && (x2 > xc+r)))
    return;                  
  end
  if (((y1 < yc-r) && (y2 < yc-r)) || ((y1 > yc+r) && (y2 > yc+r)))
    return;                  
  end
  [t1, t2, vis] = lineCircleInt(x1,y1,x2,y2,xc,yc,r);  
  if (vis == 0)
    return;  
  end
  if ((t2 < 0.0) || (t1 > 1.0))
    visible = false; 
    return;
  end
  if (t1 > 0.0) 
    cx1 = x1 + (x2-x1) * t1; 
    cy1 = y1 + (y2-y1) * t1; 
  end
  if (t2 < 1.0)  
    cx2 = x1 + (x2-x1) * t2; 
    cy2 = y1 + (y2-y1) * t2; 
  end
  visible = true;
end

function [grid] = gridKamb(x, ngrid, sigma, opts)
% GRIDKAMB  Calculates grid of density estimates from direction cosine 
% data. The grid is normalized to the contour units.                
  smooth = opts.smooth; 
  hemi = opts.hemi;
  directed = opts.directed;
  proj = opts.proj;
  [ndata,mdata] = size(x);
  if (opts.method == 1) % schmidt
    a = 0.01;                                
    zUnit = ndata*0.01;       
  else % kamb
    a = (sigma*sigma)/(ndata+sigma*sigma);    
    zUnit = sqrt(ndata*a*(1.0-a));            
  end
  if (opts.directed == true) % vectors
    alpha = 1.0-2.0*a;   
  else % axes 
    alpha = 1.0-a;
  end
  switch smooth
    case 1 % inverse area  
      f = 2.0/(1.0-alpha);                      
    case 2 % inverse area squared 
      f = 3.0/((1.0-alpha)*(1.0-alpha));
    case 3 % none  
      f = 1.0;
    otherwise % exponential  
      if (directed == true) % vectors
        f = 1.0 + ndata/(sigma*sigma);
        zUnit = sqrt(ndata*(f-1.0)/(4.0*f*f));  
      else % axes
        f = 2.0*(1.0 + ndata/(sigma*sigma));
        zUnit = sqrt(ndata*(f*0.5-1.0)/(f*f));
      end
  end
  dx = 2.0/(ngrid-1);                         
  grid = zeros(ngrid, ngrid);  
  xg = -1.0;
  for i = 1:ngrid
    yg = -1.0;
    for j = 1:ngrid
      y = sphereBProject(xg,yg,proj,hemi);
      for k = 1:ndata;
        d = dot(y,x(k,:));
        if (directed == false) 
          d = abs(d);
        end
        switch smooth 
          case 0
            grid(i,j) = grid(i,j) + exp(f*(d-1.0));
          case 1
            if (d >= alpha) 
              grid(i,j) = grid(i,j) + f*(d-alpha);
            end
          case 2
            if (d >= alpha) 
              grid(i,j) = grid(i,j) + f*(d-alpha)*(d-alpha);
            end
          case 3
            if (d >= alpha) 
              grid(i,j) = grid(i,j) + f; % 1.0
            end
        end % switch
      end % k
      yg = yg + dx;
    end % j
    xg = xg + dx;
  end % i
  if opts.mud 
    if directed then
      g = 2.0 * f / ndata;
    else
      g = f / ndata;
    end
  else
    g = 1.0/zUnit;  
  end
  grid = grid * g;
end

function [x, y, bool] = interpolate(x1, y1, z1, x2, y2, z2, z0)
% INTERPOLATE  Determine linear interpolation point between two nodes. 
  dz1 = z0-z1; 
  dz2 = z0-z2;
  if (dz1 == 0.0) 
    x = x1; 
    y = y1; 
    bool = true;
  elseif (dz2 == 0.0) 
    x = x2; 
    y = y2; 
    bool = false;
  elseif (((dz1 > 0.0) && (dz2 > 0.0)) || ((dz1 < 0.0) && (dz2 < 0.0))) 
    x = 0.0; 
    y = 0.0; 
    bool = false; 
  else
    dz = z2-z1;
    t = dz1/dz;
    x = x1 + (x2-x1) * t; 
    y = y1 + (y2-y1) * t;
    bool = true; 
  end
end

function [lines] = contourGrid(lines, x1, y1, x2, y2, grid, level)
% CONTOURGRID  Output one contour level by linear interpolation among grid nodes.
  [ng,mg] = size(grid);
  dnx = (x2-x1)/(ng-1.0); 
  dny = (y2-y1)/(mg-1.0);
  %z = level;
  gy1 = y1;
  nx = x1;
  for i = 1:ng-1
    ny = gy1;
    nxp = nx + dnx;
    for j = 1:mg-1
      nyp = ny + dny;
      z1 = grid(i,j); 
      z2 = grid(i+1,j);
      z3 = grid(i+1,j+1); 
      z4 = grid(i,j+1);
      found = 0;
      [x1,y1,bool] = interpolate(nx,ny,z1,nxp,ny,z2,level);
      if bool
        found = found+1;
      end
      [x2,y2,bool] = interpolate(nxp,ny,z2,nxp,nyp,z3,level);
      if bool
        found = found+2;
      end
      [x3,y3,bool] = interpolate(nxp,nyp,z3,nx,nyp,z4,level);
      if bool
        found = found+4;
      end
      [x4,y4,bool] = interpolate(nx,nyp,z4,nx,ny,z1,level);
      if bool
        found = found+8;
      end
      switch (found) 
        case  3 
          lines = cLineOut(lines,x1,y1,x2,y2); 
        case  5
          lines = cLineOut(lines,x1,y1,x3,y3); 
        case  9
          lines = cLineOut(lines,x1,y1,x4,y4); 
        case  6 
          lines = cLineOut(lines,x2,y2,x3,y3); 
        case 10 
          lines = cLineOut(lines,x2,y2,x4,y4); 
        case 12 
          lines = cLineOut(lines,x3,y3,x4,y4); 
        case 15
          d1 = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2)); 
          d2 = sqrt((x2-x3)*(x2-x3) + (y2-y3)*(y2-y3));
          d3 = sqrt((x3-x4)*(x3-x4) + (y3-y4)*(y3-y4)); 
          d4 = sqrt((x4-x1)*(x4-x1) + (y4-y1)*(y4-y1));
          if ((d1+d3) < (d2+d4)) 
            lines = cLineOut(lines,x1,y1,x2,y2); 
            lines = cLineOut(lines,x3,y3,x4,y4);      
          else 
            lines = cLineOut(lines,x2,y2,x3,y3); 
            lines = cLineOut(lines,x1,y1,x4,y4);
          end
      end % switch
      ny = nyp;
    end % j 
    nx = nxp;
  end % i
end

function [lines] = lineOut(lines, x1, y1, x2, y2) 
% LINEOUT  Output a line segment 
  lines = [lines; [x1,y1,x2,y2]];
end

function [lines] = cLineOut(lines, x1, y1, x2, y2) 
% VLINEOUT  Output a line segment clipped to current projection. 
  [cx1, cy1, cx2, cy2, visible] = clipLineCircle(0.0, 0.0, 1.0, x1, y1, x2, y2); 
  if (visible)
    lines = [lines; [cx1,cy1,cx2,cy2]];
  else
    lines = lines;
  end
end

function [grid, lines] = contour(dc, ngrid, sigma, nlevels, cint, cmin, cmax, opts) 
% CONTOUR  Grids data and outputs contours. 
  grid = gridKamb(dc, ngrid, sigma, opts);  
  zmin = 0.0;
  zmax = max(max(grid));
  x1 = -1.0; 
  y1 = -1.0;
  x2 = 1.0; 
  y2 = 1.0;
  lines = zeros(0,4);
  if (opts.cinterval == true) % set interval
    level = cmin;
    while (level < cmax+1e-9)
      lines = contourGrid(lines, x1, y1, x2, y2, grid, level);
      level = level + cint;
    end
  else % spaced levels
    zinc = (zmax-zmin)/nlevels;
    level = zmin;
    for i = 1:nlevels-1
      level = level + zinc;
      lines = contourGrid(lines, x1, y1, x2, y2, grid, level);
    end
  end
end

function [lines] = drawCircle(lines, x, y, radius, n)
% DRAWCIRCLE  Output a circle, adopted from Rodgers and Adams, 1976, p. 216. 
  ainc = 2.0 * pi/n;
  c1 = cos(ainc); 
  s1 = sin(ainc);
  x1 = x + radius; 
  y1 = y;
  for i = 0:n
    x2 = x + (x1-x)*c1 - (y1-y)*s1; 
    y2 = y + (x1-x)*s1 + (y1-y)*c1;
    lines = lineOut(lines, x1,y1,x2,y2);
    x1 = x2; 
    y1 = y2;
  end
end

function [frame] = drawFrame(frameopt)
% DRAWFRAME  Output projection frame. 
  frame = zeros(0,4);
  if (frameopt == 0)
    ts = 0.05;
    frame = lineOut(frame, 1.0, 0.0, 1.0-ts, 0.0);
    frame = lineOut(frame, -1.0, 0.0, -1.0+ts, 0.0);
    frame = lineOut(frame, 0.0, 1.0, 0.0, 1.0-ts);
    frame = lineOut(frame, 0.0,-1.0, 0.0, -1.0+ts);
  end
  if ((frameopt == 0) || (frameopt == 1))
    frame = drawCircle(frame, 0.0, 0.0, 1.0, 360);
  end 
end

function [bool] = hasOption(options, option)
  a = strfind(options, option);
  bool = (size(a) > 0);
end

function [opts] = getOptions(options)
  if hasOption(options, 'dec')
    opts.thetafmt = 1;
  elseif hasOption(options, 'str')
    opts.thetafmt = 2;
  elseif hasOption(options, 'dir')
    opts.thetafmt = 3;
  else
    opts.thetafmt = 0;
  end
  if hasOption(options, 'lat')
    opts.phifmt = 1;
  elseif hasOption(options, 'inc')
    opts.phifmt = 2;
  elseif hasOption(options, 'dip')
    opts.phifmt = 3;
  elseif hasOption(options, 'nad')
    opts.phifmt = 4;
  else
    opts.phifmt = 0;
  end  
  if hasOption(options, 'rad') % radians
    opts.angfmt = 0;
  elseif hasOption(options, 'grd') % gradians
    opts.angfmt = 2;
  else % degrees
    opts.angfmt = 1;
  end  
  if hasOption(options, 'mud')
    opts.mud = true;
  else
    opts.mud = false;
  end
  if hasOption(options, 'sch') % schmidt
    opts.method = 1;
  elseif hasOption(options, 'ncon') % none
    opts.method = 2;
  else % kamb
    opts.method = 0;
  end
  if hasOption(options, 'vec') % vectors
    opts.directed = true;
  else % axes 
    opts.directed = false;
  end
  if hasOption(options, 'sma') % inverse area  
    opts.smooth = 1;
  elseif hasOption(options, 'sms') % inverse area squared 
    opts.smooth = 2;
  elseif hasOption(options, 'nsm') % none  
    opts.smooth = 3;   
  else % exponential  
    opts.smooth = 0;
  end
  if hasOption(options, 'upper')
    opts.hemi = 1;
  else
    opts.hemi = 0;
  end
  if hasOption(options, 'stereo')
    opts.proj = 1;
  else
    opts.proj = 0;
  end
  if hasOption(options, 'cint')
    opts.cinterval = true;
  else
    opts.cinterval = false;
  end
  if hasOption(options, 'ntc')
    opts.frame = 1;
  elseif hasOption(options, 'nfr')
    opts.frame = 2;
  else
    opts.frame = 0;
  end
  if hasOption(options, 'ngd')
    opts.grid = 1;
  else
    opts.grid = 0;
  end
  if hasOption(options, 'gi0')
    opts.interp = -1.0;
  elseif hasOption(options, 'gi2')
    opts.interp = 0.5;
  elseif hasOption(options, 'gi3')
    opts.interp = 1.0/3.0;
  elseif hasOption(options, 'gi4')
    opts.interp = 0.25;
  elseif hasOption(options, 'gi5')
    opts.interp = 0.2;
  elseif hasOption(options, 'gi6')
    opts.interp = 1.0/6.0;
  elseif hasOption(options, 'gi8')
    opts.interp = 0.125;
  elseif hasOption(options, 'gi10')
    opts.interp = 0.1;
  else
    opts.interp = 0.2;  
  end
end
