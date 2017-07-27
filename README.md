# Undergrad_Fekete_Research
Main.pdf contains a report on the analysis of B-Spline interpolation using Fekete points.

Main.tex contains the LaTex code used to generate the report.

Each .xlsx contains sets of pre-generated Fekete points.
"XKnotFeketeMeshes" refers to Fekete points for B-Splines of Order X (N=X).
Each row is a set of Fekete points of Y internal knot intervals (K=Y). For each file, the first row designates K=4 and the final row designates K=256.

Each .png is a graphics file used within the report.

Each .m is Matlab code used to generate results used within the report.
A description for the function of each Matlab file is located at the beginning of the code as a comment.
For feketePoints.m, there are 2 options described in the comments. Option 1: Generate Fekete points. Option 2: Pull a set of Fekete points from the "XKnotFeketeMeshes" Excel files.
