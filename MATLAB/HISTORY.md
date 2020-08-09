% spherecontour.m History  

2.0.0 (2014-09-23)
* Initial public release.

2.0.1 (2018-04-01) 
* Edited description.
* Corrected to set zmin to zero in contour function. 

2.1.0 (2020-08-08)
* Added "axis('equal')" to test.m to keep figure square.
* Updated url and copyright. 
* Moved description to be self-documenting. 
* Added check in "test.m" for "_sd" in file name to indicate strike and dip.
* Fixed errors in documentation. 
* Added multiples of uniform density (MUD) option.
* Changed function parameters to: data,options,nlevels,ngrid,cint,sigma. 
* Parameters cmin and cmax are no longer used, contours start at cint and continue to maximum density.
* Degrees are now default units.
* Option 'int' changed to 'cint'.
* Option 'ste' changed to 'stereo'.
* Option 'ncn' changed to 'ncon'.
* Option 'up' changed to 'upper'.

