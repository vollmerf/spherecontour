# Errata

### SPHERCON.C

line 478  
grid[i][j] = (grid[i][j]-0.5)*f;
  
should read  
grid[i][j] = grid[i][j]*f;

### MODKAMB.PAS

line 385  
grid[i,j] := (grid[i,j]-0.5)*u;

should read  
grid[i,j] := grid[i,j]*u;

