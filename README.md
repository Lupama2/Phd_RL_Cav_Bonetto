# Phd_RL_Cav_Bonetto
Cavitation model developed by Bonetto et. al.


Below I summarize the functioning of the Puente code. The code is written in C++ using Borland Builder 6. It has six main files:
- ecdiff.cpp: This file defines the time-dependent differential equation.
- FUNCTIONS.cpp: This file defines basic functions.
- PSULSO.cpp: This file defines the model parameters.
- Rkqc.cpp: This file defines the numerical method called rkqc.
- rk4.cpp: This file defines the RK4 numerical method.
- sbsl.cpp: This is the main file that computes the dynamics of the cavitation bubble.
In addition, there are other files related to Borland and finally .dat files that contain the results/output of the evolution.

There is also the file PulWild.brp, from which the Borland interface is launched and which allows us to compile and run the code.

The directory containing the code must not include a hyphen (-) in its name!!

To compile the code from the command line, use:
bcb -B PulWid.bpr

To run it, execute the generated .exe file.

To compile individual files, use:
bcc32 -c PSULSOAR.cpp
