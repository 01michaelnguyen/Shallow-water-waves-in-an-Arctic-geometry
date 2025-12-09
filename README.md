# Shallow-water-waves-in-a-two-and-N-basin-geometry
This code visualises (linear) Shallow-water waves in time going around a two-basin geometry reminiscent to a simple Arctic Ocean model using MATLAB. 

This code also permits waves travelling to N-basin geometries also. Below are some plots from a two-basin and three-basin geometry moving in time t=0, T/3, 2T/3, where T is the period of the wave (the time taken for the wave to make one cycle to its original position around the basin). Here are some various modes travelling around different types of geometries.

Frequency and basin type| t=0 | t=T/3 | t=2T/3 |
| :-------------: | :-------------: | :-------------: |:-------------:|
|Two-basin, $\omega=0.2154$|<img src="images/eigenvalue_0.2154_frame1_streamfunction.png" width=100% height=100%>  |<img src="images/eigenvalue_0.2154_frame2_streamfunction.png" width=100% height=100%> |<img src="images/eigenvalue_0.2154_frame3_streamfunction.png" width=100% height=100%>|
|Two-basin, $\omega=0.3689$, small basins $\sigma_{1}=\sigma_{2}=-1.5$|<img src="images/eigenvalue_0.3689_frame1_streamfunction_sigma1_2_-1.5.png" width=100% height=100%>  |<img src="images/eigenvalue_0.3689_frame2_streamfunction_sigma1_2_-1.5.png" width=100% height=100%> |<img src="images/eigenvalue_0.3689_frame3_streamfunction_sigma1_2_-1.5.png" width=100% height=100%>|
|Two-basin, $\omega=0.2171$, asymmetric $\sigma_{1}=-0.25,\sigma_{2}=-0.5$|<img src="images/eigenvalue_0.2172_frame1_streamfunction_s1_-0.25_s2_-0.5.png" width=100% height=100%>  |<img src="images/eigenvalue_0.2172_frame2_streamfunction_s1_-0.25_s2_-0.5.png" width=100% height=100%> |<img src="images/eigenvalue_0.2172_frame3_streamfunction_s1_-0.25_s2_-0.5.png" width=100% height=100%>|
|Three-basin (joined), $\omega=0.2159$|<img src="images/eigenvalue_0.2159_n_3_frame1_streamfunction.png" width=100% height=100%>  |<img src="images/eigenvalue_0.2159_n_3_frame2_streamfunction.png" width=100% height=100%> |<img src="images/eigenvalue_0.2159_n_3_frame3_streamfunction.png" width=100% height=100%>|
|Three-basin (aligned), $\omega=0.2179$*|<img src="images/eigenvalue_0.2179_n_3_collinear_frame1_streamfunction-1.png" width=100% height=100%>  |<img src="images/eigenvalue_0.2179_n_3_collinear_frame2_streamfunction-1.png" width=100% height=100%> |<img src="images/eigenvalue_0.2179_n_3_collinear_frame3_streamfunction-1.png" width=100% height=100%>|

*not yet included in the code; will be finalised later.

# Using the code
Once the repository is downloaded, open the file "plot_contours_eigenvalues_twobasin_final.m". This plots the waves going around two- and N-basin geometries (using the open source cmocean colourmap to show positive and negative vorticity waves).

There are many parameters involved in plotting; those to calculate the eigenvalues to depict the (infinitely many) different wave modes depending on the initial guess of $\omega$. There are also the geometric parameters of the basins that can be changed which also permit different waves. We can also change the number of basins present in the geometry (in this version, all basins meet at a single point, see the table above).

## Initialisation
-Change "plottype" to change the output shown depicting how waves move in time. plottype=0 is the default and shows the streamfunction (contours of the streamline velocity paths) of the waves.
-Change "plot_modes" to change the mode type of the waves. plot_modes=0 is the default, non-trivial waves; but plot_modes=1 applies for symmetric basins where we have another, trivial mode that can also be plotted. plot_modes=2 is a special plot that superposes both the non-trivial and trivial mode together to create a "transferring" mode, also only present in symmetric basins.
-Change "plotwhat" to 0 to plot either a single wavemode (automatically set to the largest eigenvalue found from the initial guess of $\omega$) going around the basin in time, or 1 to plot all the possible wavemodes found from calculating all the eigenvalues the guess range of $\omega$.

## Calculating eigenvalues
Change omega_min

