Orbital animations
-------------------

The animations included here illustrate the nature of Earth's elliptical orbit about the sun, and the changing month lengths over the past 150 kyr related to variations in Earth's orbital elements.  The animations were created in R, and can be reproduced using the included scripts (adjusting for the local environment).  The R scripts create a set of .pngs, one per time step (or frame), and these can be converted to animated .gifs using ImageMagick.

The folders include: 
	
	/data			
		/ input     ! orbital elements and month lengths for the past 150 kyr, and a text file of 
		! months and days in a 365-day year

		Six subfolders, containing individual .pngs:

		/ orb_0ka
		/ orb_insol65n_anm
		/ orb_monangle
		/ orb_monlen
		/ orb_position
		/ orb_solsitice_anm

	/files			! text files with ImageMagick 7 command lines and file lists
	/gifs			! animated .gifs 
	/ppt			! a folder containing the Bartlein and Shafer (2019) AGU PowerPoint.
	/R    			! R scripts 

(Note that pathnames in the R scripts should be adjusted to reflect the local environment.)

The animations include:

- `orb_0ka.gif` The present-day elliptical orbit, showing the position of the Earth at 5-degree intervals during the year. 
- `orb_insol65n_anm`  Insolation long-term mean differences (from 0 ka, also known as "anomalies").
- `orb_position.gif` Earth's orbit over at 1 kyr intervals over the past 150 kyr.
- `orb_monlen.gif`  An animated version of Fig. 1 in Bartlein and Shafer (2019, *GMD*), showing the orbit, and month-length "anomalies" (long-term mean differences between 0 ka (1950 CE) and each time.
- `orb_solstice_anm`  An animated version of Fig. 2 in Bartlein and Shafer (2019), showing the orbit, and differences between the middle day of each month and the June solstice.  
- `orb_monangle.gif`  As alternative illustration of the month-lenth differences, where January 1 (as opposed to the March equinox is the reference point.  This animation illustrates the conservation of the "angular" or celestial definitions of the months.

For the sake of illustration, eccentricity has been exaggerated by ten for plotting in each animation, but the geometric positions, and the orbital speed in `orb_0ka.gif`, have been calculated using the true eccentricity values. 


Citation:  If you use the animations, please cite:  Bartlein, P.J. and S.L. Shafer, 2019, Paleo calendar effects on radiation, atmospheric circulation, and surface temperature, moisture, and energy-balance variables can produce interpretable but spurious large-scale patterns and trends in analyses of paleoclimatic simulations. PP31A-08, AGU 2019 Fall Meeting.  [[https://agu.confex.com/agu/fm19/meetingapp.cgi/Paper/525140]](https://agu.confex.com/agu/fm19/meetingapp.cgi/Paper/525140), as will as the GMD paper.

