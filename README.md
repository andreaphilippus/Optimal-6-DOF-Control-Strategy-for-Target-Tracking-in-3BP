# Optimal 6 Degree-of-Freedom Control Strategy for Target Tracking in Three Body Problem

https://github.com/andreaphilippus/Optimal-6-DOF-Control-Strategy-for-Target-Tracking-in-3BP/

Begun in Summer 2022 and continued in Fall 2022 as an individual study project (AAE 59700ZS) at Purdue University (West Lafayette); supervised by [Professor Kenshiro Oguri](koguri@purdue.edu).

## Problem Description

An active spacecraft is to approach and dock to a passive spacecraft in a periodic orbit around Earth-Moon L2 Lagrange point with optimal amount of fuel.

## Dependencies

The code uses [CasADi](https://web.casadi.org/) library for optimal control calculation.
* Refer to ```code/casadi-windows-matlabR2016a-v3.5.5/```.

The code requires MATLAB version of 2021b or later to run ```draw_video_CasADi.m``` as it uses ```exportgraphics``` feature which is added on version 2021b. 
* Negligible if ```draw_video_XXX``` lines in ```main.m``` are commented out.

## Author

Minyoung Andrea-Philippus (Philip) Ra, BSAAE 22' MSAA 23' @ Purdue University - West Lafayette

[GitHub](https://github.com/andreaphilippus/)

[LinkedIn](https://www.linkedin.com/in/philipra/)

[email](mra@purdue.edu)

## Acknowledgements

This project is a continuation of the final project in a course, Applied Optimal Control and Estimation (AAE 56800), instructed by Professor Inseok Hwang in Spring 2022:

[Optimal Station-Keeping of a Satellite around the Sun-Earth L1 Lagrange Point using Model Predictive Control](https://github.com/andreaphilippus/Sp22AAE568T3FinalProject)



## Contributors

Parts of the code are written/contributed by the members of the final project team in AAE 56800:

* [Mark Batistich](https://github.com/MarkBatistich) 	- [LinkedIn](https://www.linkedin.com/in/mark-batistich/)
* [Andrew Kim](https://github.com/akimb) 		- [LinkedIn](https://www.linkedin.com/in/andrewkim101/)
* [Joseph Le](https://github.com/josephtule) 		- [LinkedIn](https://www.linkedin.com/in/joseph-le-844823170/)
