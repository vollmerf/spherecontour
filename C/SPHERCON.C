/*
*  Program  : Sphere Contour
*  File     : sphercon.c
*  Purpose  : Plotting and contouring of spherical orientation data
*  Language : C
*  Compiler : Borland Turbo C++ 1.00
*  Author   : F.W. Vollmer
*  Update   : 6/92, 7/93, 9/93
*
*  Portability Notes - Non-ANSI Functions
*  --------------------------------------
*  clrscr()          - clears screen in text mode
*  fnsplit()         - splits a DOS filename
*  fnmerge()         - merges a DOS filename
*  getch(), getche() - gets character from keyboard, with echo
*  DoneGraphics()    - contains all code to shut down graphics system
*  InitGraphics()    - contains all code to initialize graphics system
*  LineOut()         - contains all code to draw line in mm units
*/

#include <conio.h>     /* Turbo C++ Library: Console IO      */
#include <dir.h>       /* Turbo C++ Library: DOS directories */
#include <graphics.h>  /* Turbo C++ Library: BGI graphics    */
#include <ctype.h>     /* ANSI C Libraries...                */
#include <math.h>
#include <stdio.h>
#include <string.h>

char szInfo[] =
  " SPHERE CONTOUR - Spherical Orientation Data Contouring\n"
  " F.W. Vollmer, 1993, State University of New York, New Paltz, New York 12561\n"
  " ---------------------------------------------------------------------------\n"
  " Sphere Contour contours three dimensional orientation data, data that can\n"
  " be represented as points on a unit sphere. This includes unit axes, vectors,\n"
  " and poles to planes. Contouring is done using modified Kamb or Schmidt\n"
  " methods. Scatter plots and rotation are options. Output is to the screen or\n"
  " a DXF file. Data must be entered into a text file in one of these formats:\n"
  "\n"
  "   FORMAT    COMPONENTS                                     EXAMPLE\n"
  "   strike    strike, dip, dip octant                        034 36 SE\n"
  "   dip       dip, dip azimuth                               35 128\n"
  "   line      plunge, trend (inclination, declination)       68 122\n"
  "   sphere    colatitude, longitude                          120 345\n"
  "\n"
  " All angles are measured in degrees. Dip and plunge angles are measured\n"
  " downward from horizontal (the XY plane). Colatitude is measured downward\n"
  " from Z (vertical). Horizontal angles strike, dip azimuth and trend are\n"
  " measured clockwise from North (the Y axis). Longitude is measured\n"
  " counterclockwise from X. Each data point must occupy one line, with the\n"
  " measurements separated by spaces or commas.\n"
  "\n"
  " Press ENTER to continue...";

#define DTOR     0.01745329252  /* degrees to radians                */
#define RTOD     57.2957795131  /* radians to degrees                */
#define MAXDATA           1600  /* maximum number of data points     */
#define MAXGRID             50  /* maximum number of data points     */
#define MAXSTR              80  /* maximum size of user input string */
#define WIDTHSTR            48  /* format width for prompt strings   */

#define round(x) (int)floor((x)+0.5)
#define sqr(x)   ((x)*(x))

#define NONE     0
enum boolean     {FALSE, TRUE};
enum datatypes   {AXES, VECTORS};
enum devices     {BGI, DXF};
enum formats     {STRIKE, DIP, LINE, SPHERE};
enum hemispheres {LOWER, UPPER};
enum methods     {KAMB, SCHMIDT};
enum projections {EQUALAREA, STEREOGRAPHIC};
enum plots       {CONTOUR, SCATTER, BOTH};
enum symbols     {CROSS=1, TRIANGLE, SQUARE, HEXAGON};
enum smoothing   {INVAREA=1, INVAREASQR};

typedef struct {
  double  ci;                /* contour interval              */
  char    dataFile[MAXSTR];  /* data file name                */
  int     device;            /* BGI || DXF                    */
  int     dataType;          /* AXIES || VECTORS              */
  int     format;            /* STRIKE || DIP || PLUNGE       */
  int     hemi;              /* LOWER || UPPER                */
  int     method;            /* KAMB || SCHMIDT               */
  double  minimum;           /* minimum contour               */
  int     nGrid;             /* number of grid nodes          */
  double  netX;              /* x coordinate of center        */
  double  netY;              /* y coordinate of center        */
  char    outFile[MAXSTR];   /* DFX file name                 */
  int     plot;              /* CONTOUR || SCATTER || BOTH    */
  int     proj;              /* EQUALAREA || STEREOGRAPHIC    */
  double  radius;            /* radius                        */
  double  rot1;              /* first rotation                */
  int     rot1Axis;          /* first rotation axis           */
  double  rot2;              /* second rotation               */
  int     rot2Axis;          /* second rotation axis          */
  double  rot3;              /* third rotation                */
  int     rot3Axis;          /* third rotation axis           */
  double  sigma;             /* binomial sigma value          */
  int     smooth;            /* NONE || INVAREA || INVAREASQR */
  int     symbol;            /* NONE || CROSS || HEXAGON      */
  double  symSize;           /* symbol size in mm             */
  } optiontype;

/***  Global Variables  ***/

optiontype  opt;
double      data[MAXDATA][3];        /* data direction cosines */
int         nData;                   /* number of data points  */
double      grid[MAXGRID][MAXGRID];  /* contour grid           */
FILE       *dxf;                     /* DXF text file          */

double dev_xRatio;        /* X device/mm ratio, negative for right origin */
double dev_yRatio;        /* Y device/mm ratio, negative for top origin   */
double dev_xOrigin;       /* left X device coordinate                     */
double dev_yOrigin;       /* bottom Y device coordinate                   */
char   dev_path[MAXPATH]; /* path to graphics device driver               */

/*** User Input ***/

void GetInt(char *pmt, int *i)
  {
  char buf[MAXSTR];
  int  j;
  printf(" %-*s %12d: ",WIDTHSTR,pmt,*i);
  fgets(buf,MAXSTR,stdin);
  if (sscanf(buf,"%d",&j) == 1) *i = j;
  }

void GetDbl(char *pmt, double *x)
  {
  char   buf[MAXSTR];
  double y;
  printf(" %-*s %12g: ",WIDTHSTR,pmt,*x);
  fgets(buf,MAXSTR,stdin);
  if (sscanf(buf,"%lf",&y) == 1) *x = y;
  }

void GetStr(char *pmt, char *s)
  {
  char buf[MAXSTR];
  char t[MAXSTR];
  printf(" %-*s %12s: ",WIDTHSTR,pmt,s);
  fgets(buf,MAXSTR,stdin);
  if (sscanf(buf,"%s",t) == 1) strcpy(s,t);
  }

void ErrorMsg(char *s, int i)
  {
  if (i > 0) fprintf(stderr," %s %d. Press ENTER to continue...",s,i);
  else fprintf(stderr," %s. Press ENTER to continue...",s);
  if (getchar()) ; /* wait */
  }

/*** Rotations ***/

void maxisrot3(int axis, double theta, double t[3][3])
  {
  int a1,a2,i,j;
  double c,s;
  for (i=0; i<3; i++) for (j=0; j<3; j++) t[i][j] = 0.0;
  t[axis][axis] = 1.0;
  a1 = (axis+1) % 3;	   /* x=012 y=120 z=201 */
  a2 = (a1+1) % 3;
  c = cos(theta);
  s = sin(theta);
  t[a1][a1] = c;
  t[a2][a2] = c;
  t[a1][a2] = -s;
  t[a2][a1] = s;
  }

void mmult3(double x[3][3], double y[3][3], double z[3][3])
  {
  int i,j,k;
  for (i=0; i<3; i++) for (j=0; j<3; j++) {
    z[i][j] = 0.0;
    for (k=0; k<3; k++) z[i][j] += x[i][k] * y[k][j];
    }
  }

void GetRotMat(double r[3][3])
  {
  int    i,j;
  double s[3][3],t[3][3];
  for (i=0; i<3; i++) {for (j=0; j<3; j++) r[i][j] = 0.0; r[i][i] = 1.0;}
  if (opt.rot1Axis) maxisrot3(opt.rot1Axis-1,opt.rot1*DTOR,r);
  if (opt.rot2Axis) {
    for (i=0; i<3; i++) for (j=0; j<3; j++) t[i][j] = r[i][j];
    maxisrot3(opt.rot2Axis-1,opt.rot2*DTOR,s);
    mmult3(s,t,r);
    }
  if (opt.rot3Axis) {
    for (i=0; i<3; i++) for (j=0; j<3; j++) t[i][j] = r[i][j];
    maxisrot3(opt.rot3Axis-1,opt.rot3*DTOR,s);
    mmult3(s,t,r);
    }
  }

/*** Conversions ***/

int OctantVal(char *s, double *r)
  {
  int i;
  for (i=0; s[i] != '\0'; i++) toupper(s[i]);
  if (strcmp(s,"N") == 0)       *r =   0.0;
  else if (strcmp(s,"NE") == 0) *r =  45.0;
  else if (strcmp(s,"E" ) == 0) *r =  90.0;
  else if (strcmp(s,"SE") == 0) *r = 135.0;
  else if (strcmp(s,"S" ) == 0) *r = 180.0;
  else if (strcmp(s,"SW") == 0) *r = 225.0;
  else if (strcmp(s,"W" ) == 0) *r = 270.0;
  else if (strcmp(s,"NW") == 0) *r = 315.0;
  else return FALSE;
  return TRUE;
  }

void PTToDC(double p, double t, double dc[3])
  /* Plunge, trend in degrees to XYZ dir cos. */
  {
  double cp;
  p = p*DTOR; t = t*DTOR;
  cp = cos(p);
  dc[0] = cp*sin(t);
  dc[1] = cp*cos(t);
  dc[2] = -sin(p);
  }

int SphereProject(double dc[3], double *x, double *y,
                  int proj, int hemi, int datatype)
  {
  int    i;
  double f,t[3];
  for (i=0; i<3; i++) t[i] = dc[i];
  if (hemi == LOWER) t[2] = -t[2];
  if (datatype == AXES && t[2] < 0.0) for (i=0; i<3; i++) t[i] = -t[i];
  if (t[2] < 0.0) return FALSE;
  if (proj == STEREOGRAPHIC) f = 1.0/(1.0+t[2]);
  else f = 1.0/sqrt(1.0+t[2]);
  *x = f*t[0]; *y = f*t[1];
  return TRUE;
  }

void SphereBProject(double x, double y, double dc[3], int proj, int hemi)
  {
  double r2,f;
  r2 = (x*x)+(y*y);
  if (proj == STEREOGRAPHIC) {
    dc[2] = (1.0-r2)/(1.0+r2);
    f = 1.0+dc[2];
    }
  else {
    f = sqrt(fabs(2.0-r2));
    dc[2] = 1.0-r2;
    }
  dc[0] = f*x;
  dc[1] = f*y;
  if (hemi == LOWER) dc[2] = -dc[2];
  }

/*** System Dependent Graphics ***/

int InitGraphics(void)
  {
  int  grmode=0,grdriver=DETECT;
  if (opt.device == DXF) {
    if ((dxf = fopen(opt.outFile,"wt")) == NULL) {
      ErrorMsg("Can not open DXF file",-1);
      return FALSE;
      }
    printf(" Plotting to %s...",opt.outFile);
    fprintf(dxf,"  0\nSECTION\n  2\nENTITIES\n");
    }
  else {
    initgraph(&grdriver,&grmode,dev_path); /* driver in exe directory */
    if (graphresult() != grOk) {
      ErrorMsg("Graphics error, required BGI driver file not found",-1);
      return FALSE;
      }
    setgraphmode(getmaxmode());
    setfillstyle(SOLID_FILL,getmaxcolor());
    bar(0,0,getmaxx(),getmaxy());
    setcolor(0);
    dev_xRatio = 2.5 * (getmaxx()+1.0)/640.0;   /* use 14" VGA as model, */
    dev_yRatio = -2.5 * (getmaxy()+1.0)/480.0;  /* it has 2.5 pixels/mm  */
    dev_xOrigin = 0.0;                          /* left                  */
    dev_yOrigin = getmaxy();                    /* bottom                */
    }
  return TRUE;
  }

void DoneGraphics(void)
  {
  if (opt.device == DXF) {
    fprintf(dxf,"  0\nENDSEC\n  0\nEOF\n");
    fclose(dxf);
    }
  else {
    if (getchar());  /* wait */
    closegraph();
    }
  }

void LineOut(double x1, double y1, double x2, double y2,char *layer)
  {
  if (opt.device == DXF)
    fprintf(dxf,"  0\nLINE\n  8\n%s\n 10\n%g\n 20\n%g\n 11\n%g\n 21\n%g\n",
            layer,x1,y1,x2,y2);
  else { /* SCREEN */
    x1 = x1*dev_xRatio+dev_xOrigin;
    y1 = y1*dev_yRatio+dev_yOrigin;
    x2 = x2*dev_xRatio+dev_xOrigin;
    y2 = y2*dev_yRatio+dev_yOrigin;
    line(round(x1),round(y1),round(x2),round(y2));
    }
  }

/*** Non-System Dependent Graphics ***/

int LineCircleInt(double x1, double y1, double x2, double y2,
                  double xc, double yc, double r, double *t1, double *t2)
  /* Adopted from Rankin 1989,p.220. */
  {
  double t,a,b,c,d,disc,dxc,dyc,dx,dy;
  *t1 = *t2 = -1.0;
  dx = x2-x1; dy = y2-y1;
  dxc = x1-xc; dyc = y1-yc;
  a = dx*dxc + dy*dyc;
  b = dx*dx + dy*dy;
  c = dxc*dxc + dyc*dyc - r*r;
  disc = a*a - b*c;
  if (disc > 0.0 && fabs(b) > 1e-6) {
    d = sqrt(disc);
    *t1 = (-a + d)/b;
    *t2 = (-a - d)/b;
    if (*t1 > *t2) {t = *t1; *t1 = *t2; *t2 = t;}
    return TRUE;
    }
  return FALSE;
  }

int ClipLineCircle(double xc, double yc, double r,
                   double *x1, double *y1, double *x2, double *y2)
  {
  double x0,y0,t1,t2;
  if ((*x1 < xc-r && *x2 < xc-r) || (*x1 > xc+r && *x2 > xc+r) ||
      (*y1 < yc-r && *y2 < yc-r) || (*y1 > yc+r && *y2 > yc+r)) return FALSE;
  if (!LineCircleInt(*x1,*y1,*x2,*y2,xc,yc,r,&t1,&t2)) return FALSE;
  if (t2 < 0.0 || t1 > 1.0) return FALSE;
  x0 = *x1; y0 = *y1;
  if (t1 > 0.0) {
    *x1 = x0 + (*x2-x0) * t1;
    *y1 = y0 + (*y2-y0) * t1;
    }
  if (t2 < 1.0) {
    *x2 = x0 + (*x2-x0) * t2;
    *y2 = y0 + (*y2-y0) * t2;
    }
  return TRUE;
  }

void DrawCircle(double x, double y, double radius, int n,  char* layer)
  /* Adopted from Rodgers and Adams, 1976, p. 216. */
  {
  double ainc,c1,s1,x1,x2,y1,y2;
  int    i;
  ainc = 2.0*M_PI/n;
  c1 = cos(ainc); s1 = sin(ainc);
  x1 = x + radius; y1 = y;
  for (i=0; i<n; i++) {
    x2 = x + (x1-x)*c1 - (y1-y)*s1;
    y2 = y + (x1-x)*s1 + (y1-y)*c1;
    LineOut(x1,y1,x2,y2,layer);
    x1 = x2; y1 = y2;
    }
  }

void CLineOut(double x1, double y1, double x2, double y2, char *layer)
  {
  if (ClipLineCircle(opt.netX,opt.netY,opt.radius,&x1,&y1,&x2,&y2))
    LineOut(x1,y1,x2,y2,layer);
  }

void DrawSymbol(double x, double y, int symbol, double size, char *layer)
  {
  double w,h,l;
  switch (symbol) {
    case CROSS:
      w = 0.5*size;
      CLineOut(x,y-w,x,y+w,layer);
      CLineOut(x-w,y,x+w,y,layer);
      break;
    case TRIANGLE:
      w = 0.5*size;
      h = w*0.5*sqrt(3.0);
      CLineOut(x-w,y-h,x+w,y-h,layer);
      CLineOut(x+w,y-h,x,y+h,layer);
      CLineOut(x,y+h,x-w,y-h,layer);
      break;
    case SQUARE:
      h = 0.5*size;
      CLineOut(x-h,y-h,x+h,y-h,layer);
      CLineOut(x+h,y-h,x+h,y+h,layer);
      CLineOut(x+h,y+h,x-h,y+h,layer);
      CLineOut(x-h,y+h,x-h,y-h,layer);
      break;
    case HEXAGON:
      w = 0.5*size;
      h = w*2.0/sqrt(3.0);
      l = 0.5*h;
      CLineOut(x,y+h,x-w,y+l,layer);
      CLineOut(x-w,y+l,x-w,y-l,layer);
      CLineOut(x-w,y-l,x,y-h,layer);
      CLineOut(x,y-h,x+w,y-l,layer);
      CLineOut(x+w,y-l,x+w,y+l,layer);
      CLineOut(x+w,y+l,x,y+h,layer);
      break;
    }
  }

void DrawNetFrame(char *layer)
  {
  double tickSize = 3.0,x,y;
  DrawCircle(opt.netX,opt.netY,opt.radius,100,layer);
  x = opt.netX+opt.radius; y = opt.netY;
  LineOut(x,y,x+tickSize,y,layer);
  x = opt.netX-opt.radius;
  LineOut(x,y,x-tickSize,y,layer);
  x = opt.netX; y = opt.netY+opt.radius;
  LineOut(x,y,x,y+tickSize,layer);
  y = opt.netY-opt.radius;
  LineOut(x,y,x,y-tickSize,layer);
  }

/*** Gridding ***/

void GridModKamb(double x[][3], int nData, double sigma,
                 int datatype, int projection, int hemisphere,
                 double grid[MAXGRID][MAXGRID], int nGrid,
                 double *zMin, double *zMax)
  {
  double y[3],a,countCos,d,dx,f,xg,yg,zUnit;
  int    i,j,k;
  if (opt.method == SCHMIDT) {
    a = 0.01;
    zUnit = nData*0.01;
    }
  else {
    a = (sigma*sigma)/(nData+sigma*sigma);
    zUnit = sqrt(nData*a*(1.0-a));
    }
  if (datatype == VECTORS) countCos = 1.0-2.0*a;
  else countCos = 1.0-a;
  if (opt.smooth == INVAREA) f = 2.0/(1.0-countCos);
  else f = 3.0/sqr(1.0-countCos);
  dx = 2.0/(nGrid-1);
  for (i=0; i<nGrid; i++) for (j=0; j<nGrid; j++) grid[i][j] = 0.0;
  xg = -1.0;
  for (i=0; i<nGrid; i++) {
    yg = -1.0;
    for (j=0; j<nGrid; j++) {
      SphereBProject(xg,yg,y,projection,hemisphere);
      for (k=0; k<nData; k++) {
        d = y[0]*x[k][0]+y[1]*x[k][1]+y[2]*x[k][2];  /* dot */
        if (datatype == AXES) d = fabs(d);
        if (d >= countCos) {
          if (opt.smooth == INVAREA)         grid[i][j] += f*(d-countCos);
          else if (opt.smooth == INVAREASQR) grid[i][j] += f*sqr(d-countCos);
          else                               grid[i][j] += 1.0;
          }
        }
      yg += dx;
      }
    xg += dx;
    }
  *zMin = 1e30; *zMax = 1e-30;
  f = 1.0/zUnit;
  for (i=0; i<nGrid; i++) for (j=0; j<nGrid; j++) {
    grid[i][j] = (grid[i][j]-0.5)*f;
    if (grid[i][j] > *zMax) *zMax = grid[i][j];
    if (grid[i][j] < *zMin) *zMin = grid[i][j];
    }
  }

/*** Contouring ***/

int Interpolate(double x1, double y1, double z1, double x2, double y2, double z2,
                double *x, double *y, double *z)
  {
  double dz,dz1,dz2,t;
  dz1 = *z-z1; dz2 = *z-z2;
  if (dz1 == 0.0) {
    *x = x1; *y = y1;
    return TRUE;
    }
  if (dz2 == 0.0) {
    *x = x2; *y = y2;
    return FALSE;
    }
  if ((dz1 > 0.0 && dz2 > 0.0) || (dz1 < 0.0 && dz2 < 0.0)) return FALSE;
  dz = z2-z1;
  t = dz1/dz;
  *x = x1 + (x2-x1) * t; *y = y1 + (y2-y1) * t;
  return TRUE;
  }

void GridContour(double x1, double y1, double x2, double y2,
  double grid[MAXGRID][MAXGRID], int ng, int mg, double level, char *layer)
  {
  double d1,d2,d3,d4,dnx,dny,nx,ny,nxp,nyp;
  double gy1,x3,x4,y3,y4,z,z1,z2,z3,z4;
  int    i,j,found;
  dnx = (x2-x1)/(ng-1.0); dny = (y2-y1)/(mg-1.0);
  z = level;
  gy1 = y1;
  nx = x1;
  for (i=0; i<ng-1; i++) {
    ny = gy1;
    nxp = nx + dnx;
    for (j=0; j<mg-1; j++) {
      nyp = ny + dny;
      z1 = grid[i][j];      z2 = grid[i+1][j];
      z3 = grid[i+1][j+1];  z4 = grid[i][j+1];
      found = 0;
      if (Interpolate( nx,ny,z1,nxp,ny, z2,&x1,&y1,&z)) found += 1;
      if (Interpolate(nxp,ny,z2,nxp,nyp,z3,&x2,&y2,&z)) found += 2;
      if (Interpolate(nxp,nyp,z3,nx,nyp,z4,&x3,&y3,&z)) found += 4;
      if (Interpolate( nx,nyp,z4,nx,ny, z1,&x4,&y4,&z)) found += 8;
      switch (found) {
        case  3: CLineOut(x1,y1,x2,y2,layer); break;
        case  5: CLineOut(x1,y1,x3,y3,layer); break;
        case  9: CLineOut(x1,y1,x4,y4,layer); break;
        case  6: CLineOut(x2,y2,x3,y3,layer); break;
        case 10: CLineOut(x2,y2,x4,y4,layer); break;
        case 12: CLineOut(x3,y3,x4,y4,layer); break;
        case 15:
          d1 = sqrt(sqr(x1-x2) + sqr(y1-y2)); d2 = sqrt(sqr(x2-x3) + sqr(y2-y3));
          d3 = sqrt(sqr(x3-x4) + sqr(y3-y4)); d4 = sqrt(sqr(x4-x1) + sqr(y4-y1));
          if ((d1+d3) < (d2+d4)) {
            CLineOut(x1,y1,x2,y2,layer); CLineOut(x3,y3,x4,y4,layer);
            }
          else {
            CLineOut(x2,y2,x3,y3,layer); CLineOut(x1,y1,x4,y4,layer);
            }
        }
      ny = nyp;
      }  /* for j */
    nx = nxp;
    }  /* for i */
  }  /* GridContour */

/*** Main Procedures ***/

int LoadData(FILE *f)
  {
  char    buf[MAXSTR];
  double  ov,p,t,x[3],r[3][3];
  char    os[MAXSTR];
  nData = 0;
  GetRotMat(r);
  while (fgets(buf,MAXSTR,f)) {
    if (buf[0] == '\0') continue;
    if (opt.format == STRIKE) {
      if (sscanf(buf,"%lf %lf %s",&t,&p,os) != 3)
        if (sscanf(buf,"%lf, %lf, %s",&t,&p,os) != 3)
          if (sscanf(buf,"%lf, %lf %s",&t,&p,os) != 3) return FALSE;
      p = 90.0-p; t = t-90.0;
      if (t < 0.0) t += 360.0;
      if (!OctantVal(os,&ov)) return FALSE;
      if (fabs(ov-t) < 90.0) t += 180.0;
      }
    else {
      if (sscanf(buf,"%lf %lf",&p,&t) != 2)
        if (sscanf(buf,"%lf, %lf",&p,&t) != 2) return FALSE;
      if (opt.format == DIP) {p = 90.0-p; t += 180.0;}
      if (opt.format == SPHERE) {p -= 90.0; t = 90.0-t;}
      }
    PTToDC(p,t,x);
    data[nData][0] = x[0] * r[0][0] + x[1] * r[0][1] + x[2] * r[0][2];
    data[nData][1] = x[0] * r[1][0] + x[1] * r[1][1] + x[2] * r[1][2];
    data[nData][2] = x[0] * r[2][0] + x[1] * r[2][1] + x[2] * r[2][2];
    nData++;
    }
  return TRUE;
  }

int InputData(void) {
  FILE   *f;
  int     change,ok;
  do {
    nData = 0;
    ok = TRUE;
    clrscr();
    printf(" SPHERE CONTOUR\n");
    printf(" --------------\n");
    GetStr("Enter data file name, or X to exit program",opt.dataFile);
    if (opt.dataFile[0] == 'x' || opt.dataFile[0] == 'X') return FALSE;
    if (!(f = fopen(opt.dataFile, "rt"))) {
      ErrorMsg("File not found",-1);
      return TRUE;
      }
    GetInt("Data format (0=strike, 1=dip, 2=line, 3=sphere)",&opt.format);
    GetInt("Data type (0=axes, 1=vectors)",&opt.dataType);
    change = FALSE;
    GetInt("Enter 1 to change projection details",&change);
    if (change) {
      GetInt("Projection (0=equal area, 1=stereographic)",&opt.proj);
      GetInt("Hemisphere (0=lower, 1=upper)",&opt.hemi);
      GetDbl("X coordinate of center in millimeters",&opt.netX);
      GetDbl("Y coordinate of center in millimeters",&opt.netY);
      GetDbl("Radius in millimeters",&opt.radius);
      }
    change = FALSE;
    GetInt("Enter 1 to rotate data",&change);
    if (change) {
      GetInt("First rotation axis (1=X, 2=Y, 3=Z)",&opt.rot1Axis);
      GetDbl("First rotation angle in degrees",&opt.rot1);
      GetInt("Second rotation axis (1=X, 2=Y, 3=Z)",&opt.rot2Axis);
      GetDbl("Second rotation angle in degrees",&opt.rot2);
      GetInt("Third rotation axis (1=X, 2=Y, 3=Z)",&opt.rot3Axis);
      GetDbl("Third rotation angle in degrees",&opt.rot3);
      }
    GetInt("Plot type (0=contour, 1=scatter, 2=both)", &opt.plot);
    if (opt.plot != SCATTER) {
      GetInt("Contouring method (0=Kamb, 1=Schmidt)",&opt.method);
      if (opt.method == KAMB)
        GetDbl("Expected level for random data (1, 2 or 3)",&opt.sigma);
      GetInt("Smoothing (0=none, 1=area, 2=area sqr)",&opt.smooth);
      GetDbl("Minimum contour (0, 1 or 2 recommended)",&opt.minimum);
      GetDbl("Contour interval (1 or 2 recommended)",&opt.ci);
      GetInt("Number of grid nodes (10 to 50)",&opt.nGrid);
      }
    if (opt.plot != CONTOUR) {
      GetInt("Symbols (1=cross, 2=triangle, 3=square, 4=hex)",&opt.symbol);
      GetDbl("Symbol size in millimeters",&opt.symSize);
      }
    GetInt("Output (0=screen graphics, 1=DXF file)",&opt.device);
    if (opt.device == DXF)
      GetStr("DXF file name",opt.outFile);
    printf("\n");
    GetInt("Enter 0 to repeat option entry, or 1 to plot",&ok);
    } while (ok == FALSE);
  if (!LoadData(f)) {
    ErrorMsg("Format error in data file at line", nData+1);
    nData = 0;
    }
  fclose(f);
  clrscr();
  return TRUE;
  }

void PlotData(char *layer)
  {
  int    i;
  double xn,yn;
  if (opt.symbol == NONE) return;
  for (i=0; i<nData; i++) if (SphereProject(data[i], &xn, &yn,
  opt.proj, opt.hemi, opt.dataType)) {
    xn = opt.netX+xn*opt.radius;
    yn = opt.netY+yn*opt.radius;
    DrawSymbol(xn,yn,opt.symbol,opt.symSize,layer);
    }
  }

void Contour(char *layer)
  {
  double level,z1,z2,x1,y1,x2,y2;
  if (nData == 0) return;
  GridModKamb(data,nData,opt.sigma,opt.dataType,opt.proj,opt.hemi,grid,opt.nGrid,&z1,&z2);
  x1 = opt.netX-opt.radius;
  y1 = opt.netY-opt.radius;
  x2 = opt.netX+opt.radius;
  y2 = opt.netY+opt.radius;
  level = opt.minimum;
  while (level < z2+1e-6) {
    GridContour(x1,y1,x2,y2,grid,opt.nGrid,opt.nGrid,level,layer);
    level += opt.ci;
    }
  }

int main(int argc, char *argv[])
  {
  char drive[MAXDRIVE],dir[MAXDIR],file[MAXFILE],ext[MAXEXT];
  fnsplit(argv[0],drive,dir,file,ext);  /* get device path, assume it */
  fnmerge(dev_path,drive,dir,"","");    /* is in exe directory        */
  opt.ci      = 2.0;
  opt.minimum = 2.0;
  opt.nGrid   = 30;
  opt.netX    = 120.0;
  opt.netY    = 100.0;
  opt.plot    = 2;
  opt.radius  = 75.0;
  opt.sigma   = 3.0;
  opt.smooth  = 2;
  opt.symbol  = 4;
  opt.symSize = 2.0;
  if (argc == 2) strncpy(opt.dataFile,argv[1],MAXSTR-1);
  strcpy(opt.outFile, "out.dxf");
  clrscr();
  printf("%s", szInfo);
  if (getchar()) ; /* wait */
  while (InputData()) {
    if (nData > 0) {
      if (!InitGraphics()) break;
      DrawNetFrame("FRAME");
      if (opt.plot != CONTOUR) PlotData("DATA");
      if (opt.plot != SCATTER) Contour("CONTOUR");
      DoneGraphics();
      }
    }
  clrscr();
  printf("Sphere Contour terminated.\n");
  return 0;
  }
