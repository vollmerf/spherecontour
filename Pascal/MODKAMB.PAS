(*
*  File     : MODKAMB.PAS
*  Puropose : Modified Kamb method orientation data contouring
*  System   : Pascal, Turbo Pascal 5.5
*  Author   : F.W. Vollmer
*  Update   : 7/3/92
*)

uses
  Crt,                        {ClrScr screen clearing procedure}
  Graph;                      {Turbo Pascal BGI graphics unit}

const                         {device constants - set for 256x192 mm VGA}
  devXRatio  = 640/256.0;     {X device/mm ratio, negative for right origin}
  devYRatio  = -480/192.0;    {Y device/mm ratio, negative for top origin}
  devXOrigin = 0.0;           {left X device coordinate}
  devYOrigin = 480.0;         {bottom Y device coordinate}
  devPath    = '';            {path for Turbo Pascal BGI device driver}

const
  DR         = 0.01745329252; {degrees to radians}
  RD         = 57.2957795131; {radians to degrees}
  pi         = 3.14159265359; {pi}
  maxData    = 1000;          {maximum number of data points}
  maxGrid    = 50;            {maximum number of grid nodes}

type
  vector3    = array[1..3] of real;
  DataArray  = array[1..maxData] of vector3;
  GridArray  = array[1..maxGrid,1..maxGrid] of real;
  OptionsRec = record
    ci       : real;          {contour interval}
    dataFile : string;        {data file name}
    device   : integer;       {output device: 0 = BGI graphics, 1 = DXF}
    direct   : integer;       {0 = axial data, 1 = directed}
    format   : integer;       {0 = strike, 1 = dip, 2 = trend plunge}
    hemi     : integer;       {0 = lower hemi, 1 = upper hemi}
    method   : integer;       {0 = ModKamb, 1 = ModSchmidt}
    minimum  : real;          {minimum contour}
    nGrid    : integer;       {number of grid nodes}
    netX     : real;          {x coordinate of center}
    netY     : real;          {y coordinate of center}
    outFile  : string;        {DFX file name}
    plot     : integer;       {0 = contour, 1 = scatter, 2 = both}
    project  : integer;       {0 = equal area, 1 = stereographic}
    radius   : real;          {radius}
    sigma    : real;          {sigma confidence level}
    smooth   : integer;       {0 = no smoothing, 1 = area, 2 = area sq.}
    symbol   : integer;       {0 = none, 1 = cross, 2 = hexagon}
    symSize  : real;          {symbol size in mm}
  end;

var                           {global variables}
  x          : DataArray;     {data direction cosines}
  nData      : integer;       {number of data points}
  grid       : GridArray;     {contour grid}
  options    : OptionsRec;    {options}
  dxf        : text;          {DXF text file}

(***  User input procedures  ***)

procedure GetInt(var i: integer);
var
  s: string[10];
  err,j: integer;
begin
  Write(' <',i,'>: ');
  Readln(s);
  Val(s,j,err);                 {convert string to integer}
  if err = 0 then i := j;
end;

procedure GetReal(var x: real);
var
  s: string[20];
  err: integer;
  y: real;
begin
  Write(' <',x:0:4,'>: ');
  Readln(s);
  Val(s,y,err);                 {convert string to real}
  if err = 0 then x := y;
end;

procedure GetString(var s: string);
var
  t: string;
begin
  Write(' <',s,'>: ');
  Readln(t);
  if t <> '' then s := t;
end;

(***  Conversions  ***)

function OctantVal(oct: string; var r: real): boolean;
var
  i: integer;
begin
  OctantVal := true;
  for i := 1 to Length(oct) do oct[i] := UpCase(oct[i]);
  while Copy(oct,1,1) = ' ' do Delete(oct,1,1);
  if oct = 'N'       then r := 0
  else if oct = 'NE' then r := 45
  else if oct = 'E'  then r := 90
  else if oct = 'SE' then r := 135
  else if oct = 'S'  then r := 180
  else if oct = 'SW' then r := 225
  else if oct = 'W'  then r := 270
  else if oct = 'NW' then r := 315
  else OctantVal := false;
end;

procedure PTToDC(p,t: real; var dc: vector3);
  {Trend, plunge in degrees to XYZ dir cos.}
var
  ci: real;
begin
  p := p*DR; t := t*DR;
  ci := Cos(p);
  dc[1] := ci*Sin(t);
  dc[2] := ci*Cos(t);
  dc[3] := -Sin(p)
end;

function SphereProject(dc: vector3; var x,y: real;
  project,hemi,direct: integer): boolean;
var
  i: integer;
  f,theta,phi: real;
begin
  if hemi = 0 then dc[3] := -dc[3];
  if direct = 0 then
    if dc[3] < 0.0 then for i := 1 to 3 do dc[i] := -dc[i];
  if dc[3] < 0.0 then SphereProject := false
  else begin
    if project = 1 then f := 1.0/(1.0+dc[3])
    else f := 1.0/Sqrt(1.0+dc[3]);
    x := f*dc[1]; y := f*dc[2];
    SphereProject := true;
  end;
end;

procedure SphereBProject(x,y: real; var dc: vector3; project,hemi: integer);
  {Allows back projection of points off net for gridding.}
var
  th,ph,r2,f: real;
begin
  r2 := Sqr(x)+Sqr(y);
  if project = 1 then begin       {stereographic}
    dc[3] := (1.0-r2)/(1.0+r2);
    f := 1.0+dc[3];
  end else begin                  {equal area}
    f := Sqrt(2.0-r2);
    dc[3] := 1.0-r2;
  end;
  dc[1] := f*x;
  dc[2] := f*y;
  if hemi = 0 then dc[3] := -dc[3];
end;

(***  System dependent graphics routines  ***)

procedure InitGraphics;
var
  grDriver,grMode,errCode: integer;
begin
  if options.device = 1 then begin  {DXF}
    Assign(dxf,options.outFile);
    Rewrite(dxf);
    Writeln(dxf,0:3);
    Writeln(dxf,'SECTION');
    Writeln(dxf,2:3);
    Writeln(dxf,'ENTITIES');
  end else begin                    {BGI}
    grDriver := Detect;
    InitGraph(grDriver,grMode,devPath);
    errCode := GraphResult;
    if errCode <> grOk then begin
      Writeln(' Graphics error: ',GraphErrorMsg(errCode));
      Halt;
    end;
  end;
end;

procedure DoneGraphics;
begin
  if options.device = 1 then begin  {DXF}
    Writeln(dxf,0:3);
    Writeln(dxf,'ENDSEC');
    Writeln(dxf,0:3);
    Writeln(dxf,'EOF');
    Close(dxf);
  end else begin                    {BGI}
    Readln;
    CloseGraph;  ClrScr;
  end;
end;

procedure LineOut(x1,y1,x2,y2: real; layer: integer);
begin
  if options.device = 1 then begin  {DXF}
    Writeln(dxf,0:3);
    Writeln(dxf,'LINE');
    Writeln(dxf,8:3);
    Writeln(dxf,layer);
    Writeln(dxf,10:3);
    Writeln(dxf,x1:0:4);
    Writeln(dxf,20:3);
    Writeln(dxf,y1:0:4);
    Writeln(dxf,11:3);
    Writeln(dxf,x2:0:4);
    Writeln(dxf,21:3);
    Writeln(dxf,y2:0:4);
  end else begin                    {BGI}
    x1 := x1*devXRatio+devXOrigin;
    y1 := y1*devYRatio+devYOrigin;
    x2 := x2*devXRatio+devXOrigin;
    y2 := y2*devYRatio+devYOrigin;
    Line(Round(x1),Round(y1),Round(x2),Round(y2));
  end;
end;

(***  Non-system dependent graphics routines  ***)

function LineCircleInt(x1,y1,x2,y2: real; xc,yc,r: real;
  var t1,t2: real): boolean;
  {Adopted from Rankin 1989, p.220.}
var
  t,a,b,c,d,disc,dxc,dyc,dx,dy: real;
  n: integer;
begin
  LineCircleInt := false;
  t1 := -1; t2 := -1;
  dx := x2-x1; dy := y2-y1;
  dxc := x1-xc; dyc := y1-yc;
  a := dx*dxc + dy*dyc;
  b := Sqr(dx) + Sqr(dy);
  c := Sqr(dxc) + Sqr(dyc) - Sqr(r);
  disc := Sqr(a) - b*c;
  if (disc > 0.0) and (Abs(b) > 1e-6) then begin
    d := Sqrt(disc);
    t1 := (-a + d)/b;
    t2 := (-a - d)/b;
    if t1 > t2 then begin
      t := t1; t1 := t2; t2 := t;
    end;
    LineCircleInt := true;
  end;
end;

function ClipLineCircle(xc,yc,r: real; var x1,y1,x2,y2: real): boolean;
var
  x0,y0,t1,t2: real;
  n: integer;
begin
  ClipLineCircle := false;
  if (x1 < xc-r) and (x2 < xc-r) then Exit;
  if (x1 > xc+r) and (x2 > xc+r) then Exit;
  if (y1 < yc-r) and (y2 < yc-r) then Exit;
  if (y1 > yc+r) and (y2 > yc+r) then Exit;
  if not LineCircleInt(x1,y1,x2,y2,xc,yc,r,t1,t2) then Exit;
  if (t2 < 0.0) or (t1 > 1.0) then Exit;
  x0 := x1; y0 := y1;
  if t1 > 0.0 then begin
    x1 := x0 + (x2-x0) * t1;
    y1 := y0 + (y2-y0) * t1;
  end;
  if t2 < 1.0 then begin
    x2 := x0 + (x2-x0) * t2;
    y2 := y0 + (y2-y0) * t2;
  end;
  ClipLineCircle := true;
end;

procedure DrawCircle(x,y,radius: real; n: integer);
  {Adopted from Rodgers and Adams, 1976, p. 216.}
var
  ainc,c1,s1,x1,x2,y1,y2: real; i: integer;
begin
  ainc := 2.0*pi/n;
  c1 := Cos(ainc); s1 := Sin(ainc);
  x1 := x + radius; y1 := y;
  for i := 1 to n do begin
    x2 := x + (x1-x)*c1 - (y1-y)*s1;
    y2 := y + (x1-x)*s1 + (y1-y)*c1;
    LineOut(x1,y1,x2,y2,0);
    x1 := x2; y1 := y2;
  end;
end;

procedure DrawNetLine(x1,y1,x2,y2: real; layer: integer);
begin
  with options do
    if ClipLineCircle(netX,netY,radius,x1,y1,x2,y2) then
      LineOut(x1,y1,x2,y2,layer);
end;

procedure DrawSymbol(x,y: real; symbol: integer; size: real);
var
  h,k,l,s: real;
begin
  s := 0.5*size;
  if symbol = 1 then begin           {cross}
    DrawNetLine(x,y-s,x,y+s,2);
    DrawNetLine(x-s,y,x+s,y,2);
  end else if symbol = 2 then begin  {hexagon}
    h := s*Sin(60.0*DR);
    l := 0.5*s;
    DrawNetLine(x+s,y,x+l,y+h,2);
    DrawNetLine(x+l,y+h,x-l,y+h,2);
    DrawNetLine(x-l,y+h,x-s,y,2);
    DrawNetLine(x-s,y,x-l,y-h,2);
    DrawNetLine(x-l,y-h,x+l,y-h,2);
    DrawNetLine(x+l,y-h,x+s,y,2);
  end;
end;

procedure DrawNetFrame;
const
  tickSize = 3.0;
var
  x,y: real;
begin
  with options do begin
    DrawCircle(netX,netY,radius,100);
    x := netX+radius; y := netY;
    LineOut(x,y,x+tickSize,y,0);
    x := netX-radius;
    LineOut(x,y,x-tickSize,y,0);
    x := netX; y := netY+radius;
    LineOut(x,y,x,y+tickSize,0);
    y := netY-radius;
    LineOut(x,y,x,y-tickSize,0);
  end;
end;

(***  Gridding  ***)

procedure GridModKamb(var x: DataArray; nData: integer;
  sigma: real; direct,project,hemi: integer;
  var grid: GridArray; nGrid: integer; var zMin,zMax: real);
var
  y: vector3;
  a,countCos,d,dx,f,u,xg,yg,zUnit: real;
  i,j,k: integer;
begin
  if options.method = 1 then begin
    a := 0.01;
    zUnit := nData*0.01;
  end else begin
    a := Sqr(sigma)/(nData+Sqr(sigma));
    zUnit := Sqrt(nData*a*(1.0-a));
  end;
  if direct = 1 then countCos := 1.0-2.0*a
  else countCos := 1.0-a;
  if options.smooth = 1 then f := 2.0/(1.0-countCos)
  else f := 3.0/Sqr(1.0-countCos);
  dx := 2.0/(nGrid-1);
  FillChar(grid,SizeOf(grid),0);  {initialize grid to zero}
  xg := -1;
  for i := 1 to nGrid do begin
    yg := -1;
    for j := 1 to nGrid do begin
      SphereBProject(xg,yg,y,project,hemi);
      for k := 1 to nData do begin
        d := y[1]*x[k,1]+y[2]*x[k,2]+y[3]*x[k,3];  {dot product}
        if direct = 0 then d := Abs(d);
        if d >= countCos then
          if options.smooth = 2 then
            grid[i,j] := grid[i,j] + f*Sqr(d-countCos)
          else if options.smooth = 1 then
            grid[i,j] := grid[i,j] + f*(d-countCos)
          else
            grid[i,j] := grid[i,j]+1;
      end;
      yg := yg + dx;
    end;
    xg := xg + dx;
  end;
  u := 1.0/zUnit;
  zMin := 1e30;
  zMax := -1e30;
  for i := 1 to nGrid do for j := 1 to nGrid do begin
    grid[i,j] := (grid[i,j]-0.5)*u;
    if grid[i,j] < zMin then zMin := grid[i,j];
    if grid[i,j] > zMax then zMax := grid[i,j];
  end;
end;

(***  Contouring  ***)

function Interpolate(x1,y1,z1,x2,y2,z2: real; var x,y,z: real): boolean;
  {Find the interpolated point of value z between two nodes.}
var
  dz,dz1,dz2,t: real;
begin
  dz1 := z-z1; dz2 := z-z2;
  if dz1 = 0.0 then begin
    x := x1; y := y1;
    Interpolate := false;
  end else if dz2 = 0.0 then begin
    x := x2; y := y2;
    Interpolate := false;
  end else if ((dz1 > 0.0) and (dz2 > 0.0))
  or ((dz1 < 0.0) and (dz2 < 0.0)) then
    Interpolate := false
  else begin
    dz := z2-z1;
    t := dz1/dz;
    x := x1 + (x2-x1) * t; y := y1 + (y2-y1) * t;
    Interpolate := true;
  end;
end;

procedure GridContour(x1,y1,x2,y2: real; var grid: GridArray;
  ng,mg: integer; level: real);
var
  d1,d2,d3,d4,dnx,dny,nx,ny,nxp,nyp: real;
  gy1,x3,x4,y3,y4,z,z1,z2,z3,z4: real;
  i,j,found: integer;
begin
  dnx := (x2-x1)/(ng-1.0);
  dny := (y2-y1)/(mg-1.0);
  z := level;
  gy1 := y1;
  nx := x1;
  for i := 1 to ng-1 do begin
    ny := gy1;
    nxp := nx + dnx;
    for j := 1 to mg-1 do begin
      nyp := ny + dny;
      z1 := grid[i,j];
      z2 := grid[i+1,j];
      z3 := grid[i+1,j+1];
      z4 := grid[i,j+1];
      found := 0;
      if Interpolate(nx,ny,z1,nxp,ny,z2,x1,y1,z) then Inc(found,1);
      if Interpolate(nxp,ny,z2,nxp,nyp,z3,x2,y2,z) then Inc(found,2);
      if Interpolate(nxp,nyp,z3,nx,nyp,z4,x3,y3,z) then Inc(found,4);
      if Interpolate(nx,nyp,z4,nx,ny,grid[i,j],x4,y4,z) then Inc(found,8);
      case found of
        3  : DrawNetLine(x1,y1,x2,y2,1);
        5  : DrawNetLine(x1,y1,x3,y3,1);
        9  : DrawNetLine(x1,y1,x4,y4,1);
        6  : DrawNetLine(x2,y2,x3,y3,1);
        10 : DrawNetLine(x2,y2,x4,y4,1);
        12 : DrawNetLine(x3,y3,x4,y4,1);
        15 : begin
          d1 := Sqrt(Sqr(x1-x2) + Sqr(y1-y2));
          d2 := Sqrt(Sqr(x2-x3) + Sqr(y2-y3));
          d3 := Sqrt(Sqr(x3-x4) + Sqr(y3-y4));
          d4 := Sqrt(Sqr(x4-x1) + Sqr(y4-y1));
          if (d1+d3) < (d2+d4) then begin
            DrawNetLine(x1,y1,x2,y2,1);
            DrawNetLine(x3,y3,x4,y4,1);
          end else begin
            DrawNetLine(x2,y2,x3,y3,1);
            DrawNetLine(x1,y1,x4,y4,1);
          end;
        end;
      end;
      ny := nyp;
    end; {for j}
    nx := nxp;
  end; {for i}
end; {GridContour}

(***  Main procedures  ***)

procedure Info;
begin
  ClrScr;  {clear screen}
  Writeln;
  Writeln('                                 MODKAMB');
  Writeln;
  Writeln('                          SPHERICAL ORIENTATION');
  Writeln('                         DATA CONTOURING PROGRAM');
  Writeln;
  Writeln('                               F.W. Vollmer');
  Writeln('                       State University of New York');
  Writeln('                         The College at New Paltz');
  Writeln('                        New Paltz, New York 12561');
  Writeln('                       ****************************');
  Writeln;
  Writeln(' Enter the orientation data into a text file using either strike (strike,');
  Writeln(' dip, and dip octant), dip (dip and dip azimuth), or line format (plunge');
  Writeln(' and trend, or inclination and declination). Each data point should');
  Writeln(' occupy one line, with the measurements separated by spaces.');
  Writeln;
end;

procedure LoadData(var f: text);
  {Read data from file f and convert to dir cos.}
var
  ov,p,t: real;
  os: string[3];
begin
  nData := 0;
  while not EOF(f) do begin
    if options.format = 0 then begin
      Readln(f,t,p,os);
      p := 90.0-p; t := t-90.0;
      if t < 0.0 then t := t+360.0;
      if OctantVal(os,ov) then
        if Abs(ov-t) < 90.0 then t := t+180.0
    end else if options.format = 1 then begin
      Readln(f,p,t);
      p := 90.0-p; t := t+180.0;
    end else
      Readln(f,p,t);
    nData := nData+1;
    PTToDC(p,t,x[nData]);
  end;
end;

function InputData: boolean;
  {Get user input, read data from file and convert to dir cos.}
var
  f: text;            {text file}
  change,ok: integer;
begin
  repeat
    nData := 0;
    InputData := false;
    ok := 1;
    Writeln;
    Writeln(' Press ENTER for defaults, or enter the desired values.');
    Write(' Enter a SPACE to quit, or the data file name    ');
    GetString(options.dataFile);
    if (options.dataFile = ' ') or (options.dataFile = '') then Exit;
    Assign(f,options.dataFile); {assign file name}
    {$I-}
    Reset(f); {open file for reading}
    {$I+}
    InputData := true;
    if IOResult <> 0 then begin {input/output error}
      Writeln(' File not found.');
      Exit;
    end;
    Write(' Data format (0 = strike, 1 = dip, 2 = line)     ');
    GetInt(options.format);
    Write(' Data type (0 = axial data, 1 = vectoral data)   ');
    GetInt(options.direct);
    change := 0;
    Write(' Enter 1 to change projection details            ');
    GetInt(change);
    if change = 1 then begin
      Write(' Projection (0 = equal area, 1 = stereographic)  ');
      GetInt(options.project);
      Write(' Hemisphere (0 = lower, 1 = upper)               ');
      GetInt(options.hemi);
      Write(' X coordinate of center in millimeters           ');
      GetReal(options.netX);
      Write(' Y coordinate of center in millimeters           ');
      GetReal(options.netY);
      Write(' Radius in millimeters                           ');
      GetReal(options.radius);
    end;
    Write(' Plot type (0 = contour, 1 = scatter, 2 = both)  ');
    GetInt(options.plot);
    if options.plot <> 1 then begin
      Write(' Gridding method (0 = ModKamb, 1 = ModSchmidt)   ');
      GetInt(options.method);
      if options.method = 0 then begin
        Write(' Sigma confidence level (1 to 3 suggested)       ');
        GetReal(options.sigma);
      end;
      Write(' Smoothing (0 = none, 1 = area, 2 = area squared)');
      GetInt(options.smooth);
      Write(' Minimum contour (0, 1 or 2 recommended)         ');
      GetReal(options.minimum);
      Write(' Contour interval (1 or 2 recommended)           ');
      GetReal(options.ci);
      Write(' Number of grid nodes (10 to 50)                 ');
      GetInt(options.nGrid);
    end;
    if options.plot <> 0 then begin
      Write(' Symbols (1 = crosses, 2 = hexagons)             ');
      GetInt(options.symbol);
      Write(' Symbol size in millimeters                      ');
      GetReal(options.symSize);
    end;
    Write(' Output (0 = screen graphics, 1 = DXF file)      ');
    GetInt(options.device);
    if options.device = 1 then begin
      Write(' DXF file name                                   ');
      GetString(options.outFile);
    end;
    Writeln;
    Write(' Enter 0 to repeat option entry, or 1 to plot    ');
    GetInt(ok);
  until ok = 1;
  LoadData(f);
  Close(f);
end;

procedure PlotData;
var
  i : integer;
  xn,yn : real;
begin
  if options.symbol > 0 then for i := 1 to nData do begin
    with options do
      if SphereProject(x[i],xn,yn,project,hemi,direct) then begin
        xn := options.netX+xn*options.radius;
        yn := options.netY+yn*options.radius;
        DrawSymbol(xn,yn,options.symbol,options.symSize);
      end;
  end;
end;

procedure Contour;
var
  level : real;
  x1,y1,x2,y2,z1,z2 : real;
begin
  if nData = 0 then Exit;
  with options do
    GridModKamb(x,nData,sigma,direct,project,hemi,grid,nGrid,z1,z2);
  x1 := options.netX-options.radius;
  y1 := options.netY-options.radius;
  x2 := options.netX+options.radius;
  y2 := options.netY+options.radius;
  level := options.minimum;
  while level < z2+1e-6 do begin
    GridContour(x1,y1,x2,y2,grid,options.nGrid,options.nGrid,level);
    level := level + options.ci;
  end;
end;

begin
  FillChar(options,SizeOf(options),0);  {initialize to zero}
  options.ci := 2.0;                    {set nonzero defaults}
  options.dataFile := ' ';
  options.minimum := 2.0;
  options.nGrid := 30;
  options.netX := 120.0;
  options.netY := 100.0;
  options.outFile := 'OUT.DXF';
  options.plot := 2;
  options.radius := 75.0;
  options.sigma := 3.0;
  options.smooth := 2;
  options.symbol := 2;
  options.symSize := 2;
  Info;
  while InputData do if nData > 0 then begin
    InitGraphics;
    DrawNetFrame;
    if options.plot <> 0 then PlotData;
    if options.plot <> 1 then Contour;
    DoneGraphics;
  end;
end.
