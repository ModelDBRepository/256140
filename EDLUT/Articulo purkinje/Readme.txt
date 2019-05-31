BRIEF HOW TO FOR EDLUT SIMULATOR

Windows platform

Software to work with:

-	Windows XP/7/8 x64
-	visual studio 2008 x64

Optional Software to work with:

-	Intel C++ Composer XE 2013 (in order to take advantage of OpenMP parallel task options)


Two case studies:  
-	CASE A) Case A and B use leaky integrate-and-fire (LIF) neuron models except for the PC. CASE A Cerebellar model implements a PC deatiled model. The PC model is a reduction of the PC detailed computational model presented in [1-2]and consists of a conductance-based mono-compartment model with four ionic currents.
-	cASE B) Case study B uses a LIF neuron playing the role of a PC. CASE B is implemented for comparative purposes.

Case A and B are faced to a 1Hz complete head rotation movement (1sec) that needs to be compensated by an eye-movement right in the opposite direction (in counter phase) over 2500 iterations. 
During the next 2500 iterations, the visual stimulus is given simultaneously in phase. 
That is, in the same direction of the vestibular stimulation but at a greater amplitude (10% more).
The cerebellum is then induced to re-compensate the eye movement oppositely to the direction initially followed over the first 2500 iterations. 
Finally, the initial set-up is restored over the last 3500 iterations. 
1000 extra iterations are included to ensure the stabilisation of the VOR adaptation. 
STDP is enabled at PF-PC, MF-DCN, and PC-DCN synapses. PC-DCN STDP is meant to correlate (more concretely, to reverse-correlate) the activity arriving from PC with DCN output activity. PC-DCN plasticity adjusts the firing rate at DCN as the sequence of MF groups are activated (states) rather than contributing to learning consolidation

 
Once all the needed software is installed different paths have to be modified / inserted


1)	Open the project with visual studio. Select "Resease" in Solution Configuration and "x64" in Solution Platform. Open project property pages from the current project (whether CASE A OR CASE B).

a)	In Configuration  Properties\Debugging\working Directory it is mandatory to indicate the working directory we are going to use. This directory has to contain all the needed files that EDLUT needs to properly operate (neural network definition, neuron files, weight files, etc.)
 	(in our case) ..\WORKING_DIRECTORY\WINDOWS


b) 	Set "configuration Properties\C/C++\Language\OpenMP Support" to Yes (/openmp). 

2)	COMPILE THE WHOLE PROJECT

3)	RUN THE PROJECT IN VISUAL STUDIO OR IN A TERMINAL (this option is faster):
	Open a terminal and go to the following path:
	"folder path/WORKING_DIRECTOTY/WINDOWS/"

	Run the simulation:
	"..\..\CASE_A\x64\Release\CASE_A"
	"..\..\CASE_B\x64\Release\CASE_B"


4)	EDLUT generates the file vars.dat which contains all the needed variables to be checked afterwards

5)	EDLUT also generates the file ouputactivity.dat which contains the obtained activity from the neural-network

6)	By using the script programmed in  MATLAB varlog_reduced_VOR.m  we can easily visualize all the already stored variables in vars.dat 

7)	In the same way, by using the script programmed in MATLAB noout2par.m we can easily visualize the stored activity in activity.dat

...................................................................................................................................................................



BRIEF HOW TO FOR EDLUT SIMULATOR
Linux platform

Software we are working with:

-	Compiler gcc, (It can also work with clang for Mac Os X )

Once all the needed software is installed we have to follow the following steps for CASE_A and CASE_B:
1)	Open a terminal and go to the folder CASE_A or CASE_B.

2)	Give execution permision to the configure file by typing:
	"chmod u+x configure".

3)	Create the makefile disabling matlab and cuda support: 
	"./configure --without-matlab --without-cuda".

4)	Compile the proyect:
	"make"

5)	Once the source code has been compiled, the case studies can be simulated. Open a terminal and go to the following path:
	"folder path/WORKING_DIRECTOTY/UNIX/"

6) 	Run the project:
	"../../CASE_A/bin/robottest"
	"../../CASE_B/bin/robottest"


7)	EDLUT generates the file vars.dat which contains all the needed variables to be checked afterwards

8)	EDLUT also generates the file ouputactivity.dat which contains the obtained activity from the neural-network

9)	By using the script programmed in  MATLAB varlog_reduced_VOR.m  we can easily visualize all the already stored variables in vars.dat 

10)	In the same way, by using the script programmed in MATLAB noout2par.m we can easily visualize the stored activity in activity.dat
......................................................................................................................................................................


[1].	Middleton SJ, Racca C, Cunningham MO, Traub RD, Monyer H, Knopfel T, et al. High-frequency network oscillations in cerebellar cortex. Neuron. 2008;58(5):763-74.
[2].	Miyasho T, Takagi H, Suzuki H, Watanabe S, Inoue M, Kudo Y, et al. Low-threshold potassium channels and a low-threshold calcium channel regulate Ca2+ spike firing in the dendrites of cerebellar Purkinje neurons: a modeling study. Brain research. 2001;891(1-2):106-15.
