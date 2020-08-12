# SphereContour
Program for automatic contouring of spherical orientation data using a modified Kamb method to create spherical projections including Schmidt (equal-area) and stereographic (equal-angle) plots. Additional files include example data files. 

## Citation
This repostory contains a C implementation (SPHERECON) of the algorithms presented in:

__Vollmer, F.W., 1995. C program for automatic contouring of spherical 
orientation data using a modified Kamb method: Computers & Geosciences, 
v. 21, n. 1, p. 31-49.__

which should be cited by publications using this code, algorithm, or derivative 
works, to produce figures or other content. 

## License
Please see [LICENSE](./LICENSE.md) and [CITATION](CITATION.md) for license and citation information.

## Notes
This is the C code published in Vollmer (1995). It is a command line program (SPHERECON) created using Borland Turbo C++ and BGI graphics. It has not been updated to run on current C++ compilers, such as GNU GCC, and has been superseded by [__Orient__](https://vollmerf.github.io/orient/).

## Related Software
The program [__Orient__](https://vollmerf.github.io/orient/) by this author is free software that includes this functionality, and numerous additional options. It runs on Macintosh, Windows, and Linux platforms.
