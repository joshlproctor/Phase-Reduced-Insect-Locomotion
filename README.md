# Phase-Reduced-Insect-Locomotion

Introduction:

This repository contains MATLAB scripts implementing the phase-reduced model of cockroach locomotion described in [1].   In this work, we investigate the effects of proprioceptive and exteroceptive feedback on stability and maneuverability.  The model is an extension of earlier work on phase-reduced models of legged locomotion [2] and a broader set of investigations of legged cockroach locomotion [3,4,5].  

The MATLAB scripts in this repository provide an example code base for describing the phase-reduced model of [1].  Specifically, we include the `main.m` script that simulates forward the cockroach locomotion model from an initial condition.  This code provides an example that implements the lateral perturbation experiments with varying synaptic strengths.  This corresponds to the investigations described in Section 4.1 of [1].  The folder `supportingFunctions` includes all of the helper functions and input data required for `main.m`.  Note that the `main.m` script simulates forward the model from an initial condition.  The script will save the workspace of the simulation and produce a simple plot of the center of mass trajectory of the cockroach.  Users can adapt this code for different investigations including turning or changing the mass or moment of inertia along with visualizing the output in different formats.

In addition to the code, we have also supplied a small set of documentation.   The documentation includes a file that describes how to execute main.m called `Running the code`.  There is also a file called `Description of supporting functions` which has short descriptions of each of the files in the supportingFunctions folder.  We also include the numbering convention for the neural system model which is directly implemented in the code.  The numbering convention is schematically described in the file `NamingDiagram.pdf`.  This schematic has been reproduced from an earlier repository [6] which implements the phase-reduced model in [2].

Executing main.m

The `main.m` MATLAB script simulates the phase-reduced model.  Specifically, the file implements the lateral cannon perturbation with varying synaptic feedback strengths.  For more details about running the file, please see the documentation file `Running the code`.

Supporting functions

All of the supporting functions for main.m are contained within the folder `supportingFunctions`.  In the `Documentation` folder, we have included a file with a short description for each function.  

Documentation

The `Documentation` folder contains descriptions of how to execute the code, supporting functions for running main.m, and the naming descriptions for the neural system and legs.  


[1]  Proctor, J.L., Holmes, P. The effects of feedback on stability and maneuverability of a phase-reduced model for cockroach locomotion. Biological Cybernetics (in review).  

[2]  Proctor J, Kukillaya R, Holmes P (2010) A phase-reduced neuro-mechanical model for insect locomotion: feed-forward stability and proprioceptive feedback. Phil Trans R Soc A 368:5087–5104

[3]  Schmitt J, Holmes P (2000) Mechanical models for insect locomotion: Dynamics and stability in the horizontal plane – I. Theory. Biol Cybern 83(6):501–515

[4] Seipel J, Holmes P, Full R (2004) Dynamics and stability of insect locomotion: A hexapedal model for horizontal plane motion. Biol Cybern 91(2):76–90

[5]  Kukillaya R, Holmes P (2007) A hexapedal jointed-leg model for insect locomotion in the horizontal plane.
Biol Cybern 97:379–395

[6]  Electronic Physics Auxiliary Publication Service E (2009) See document no. e-chaoeh-19-005992 for parameter values and code documentation. URL: http://ftp.aip.org/epaps/chaos/E-CHAOEH-19-005992/ For more information on EPAPS, see http://www.aip.org/pubservs/epaps.html

