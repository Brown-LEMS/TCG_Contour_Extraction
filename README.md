This is the repository of TCG curve detector(in Matlab).

Usage:

To run TCG, one needs to first generate SE edge. 

Run TCG_contour_extraction_matlab/util/SE_Dollar_detector/edgesDemo_SETO_ucm.m (Don't forget to change the path in the script), you will get one .edg file, one _bry.png file and one _SE.png file.

Run TCG_main.m. .cem file will be generated. 

Corner Detector: "corner_pts" outputed from "contour_breaker_at_conner()" function will be the final corner points detected by Yuliang. 