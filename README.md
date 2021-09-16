# SingleParticleTracking
Matlab code for MSD analysis changed based on @msdanalyzer



\section{Fiji/ImageJ: Particle Tracker 2D/3D}
     
\textbf{Open and adjust Images/videos:}

1. Drag the .nd2 file into Fiji;

2. Click OK for the first pop-up window (Bio-Formats Import Options) (For most of the usual cases, don’t have to change the selections);

3. Close the second pop-up window (Original Metadata – File name), or keep it for information check, eg. exposure time, timestamp information…; 

4. Wait for the image lookup window to show up; 

5. If image shows totally black, adjust the brightness: 

Image > adjust > Brightness/Contrast…Usually click 'Auto', then 'Apply' should be fine; if not, adjust manually.\\

\noindent
\textbf{Set the right pixel/micron ratio:}

6. Analyze > Set Scale: usually 9.375 pixels/micron for our microscope setup (60X objective and a 2.5X extra magnification in the CCD camera);
Check 'Global', click 'OK'.\\

\noindent
\textbf{Particle tracking: (The most time-consuming part)}

7. Plugins > Mosaic > Particle Tracking 2D/3D;

8. Find the optimal parameters for particle detection and linking and tracking, use 'Preview Detected' to preview the detection: Usually, Radius: 3-5; Per/Abs: 0.01-0.05; Link Range: 4; Displacement: 4-7; Dynamics: Brownian; (Refer to the attached paper for more information about the tracking parameters)

9. Click 'OK', wait for All Trajectories Visual window and Results window to show up.\\

\noindent
\textbf{Collect and save trajectory information:}

10. In the All Trajectories Visual window, go through the entire video, check if most of the trajectories of moving particles are well recognized; if not, go back to Step 9, adjust the parameters and start over again;

11. If good, go to the Results window, 

12. Click Save Full Report to save all parameter settings and motion information for this video detection;
	
13. Click All Trajectories to Table, double click > Save as > XXX\_\#.csv file (for MATLAB MSD Analysis later);
	
14. Click All MSS/MSD to Table; ( At the first time, there might be a pop-up window warning about dimension unit > click OK > change Unit of length from “micron” to “um” > check Global > click OK); 

15. If Reset Results Table window pop up, click OK; > double click > Save as > MSD\_XXX\_\#.csv file.\\ 

\noindent
For detailed tutorial, please refer to: \newline
\url{https://mosaic.mpi-cbg.de/MosaicToolboxSuite/ParticleTracker.html}



\section{MATLAB: Trajectory filtration and MSD analysis}
There are two sets of code for Trajectory filtration and MSD analysis:

- PlotLogScaleMSDForUnknownDiffGeneral.m

- PlotLogScaleMSDForBrownianDiffGeneral.m

\noindent
The only difference between these two codes is the filtration condition for filter out/recognize good trajectories. In UnknownDiff, $\big<(\Delta r_i)^2\big>=4D_\textrm{general}t^\alpha$, log($\big<(\Delta r_i)^2\big>$) and log($t$) are fitted to a linear line, defining $R^2$ > 0.9 as good enough trajectory. For BrownianDiff, $\big<(\Delta r_i)^2\big>=4D_\textrm{linear}t$, $\big<(\Delta r_i)^2\big>$ and $t$ are fitted to a linear line, also with $R^2$ > 0.9 as good enough trajectory. In both methods, you can get Dgeneral, Dlinear and anomalous exponent $\alpha$ for each good trajectory that is filtered out.\\

1. Open “PlotLogScaleMSDForUnknown/BrownianDiffGeneral.m” in MATLAB;

2. Collect all raw data files (XXX\_30ms\_\#.csv) in one folder and duplicate (always keep one as backup for other analysis maybe in the future); 
	
3. Open the folder in MATLAB as Current Folder shown on the left hand side;
	
4. Change the initial settings: 

- concentration (change this according to the file name)

- frame\_interval (get this by subtracting two continuous timestamp)

- video (give video \# range)

- remove ( to remove \# that are missing)

- clip\_factor (usually set as 0.25, only fit the first one quarter of the MSD plot)

- Mini\_Trajlength (define the mini trajectory length to rule out those short trajectories)
 
5. Click Run (green button on top);
	
6. Wait to get: (Files will be automatically generated and saved under the Current Folder selected.)

- MSD data for each individual good trajectory:

\textit{concentration$\_$video$\#\_$MSDdata$\_$Traj$\#$.xlsx}; 

- All information ($D_\textrm{general}, D_\textrm{linear},~\alpha,~R^2$ etc.) for all good trajectories under this concentration:

\textit{AllDcollections\_(un)brownian\_Traj>miniTrajlenght\_concentration.xlsx};

- A log scale MSD plot for all good trajectories: 

\textit{concentration\_logScale\_MSD ALL.png)}. 



  




