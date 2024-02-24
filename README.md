# Finite-Elements-Method
This set of scripts was a reasonable way to learn the rudiments of the Finite Elements Method(FEM). I have taken some well known problems which have analytical solution and solved them computationally. I did this because 
I wanted to make sure that my results are consistent and also acquire confidence with my skills. Thus, in every script the analytical solution (or result) is also provided so as the user can make their own tests. 

I have implemented first order scalar elements to solve the following problems: 

1) Cylindrical capacitor.
2) Parallel plate capacitor
3) Circular and Rectangular metallic waveguides
4) 2D-Scattering problems.
   
A better approximation can be realized by increasing the nubmer of refinements in the geometry. 

For instance, **Parallel_Capacitor.m** implements three different solvers and produces figures like the following 

![Parallel_1](https://github.com/Filippos-Dakis/Finite-Elements-Method/assets/114699564/d4208fac-2680-4ddc-a2db-db323b6271c4)

**Circular_Waveguide.m** script solves the EM eigenvalue problem in a circular waveguide

![Waveguide_1](https://github.com/Filippos-Dakis/Finite-Elements-Method/assets/114699564/3d31e5f1-e4cc-449e-9588-530e0523c46d)


while the **Scattering_Problem_TM_polarization.m** solves the scattering problem of a TM polarized plane wave incident on a infinite cylinder

![Scattering_1](https://github.com/Filippos-Dakis/Finite-Elements-Method/assets/114699564/c6226370-b2a3-4a56-8136-815f4821e765)



A last remark is that the scripts are written in user friendly way and I strongly encourage the user to go through the code and read it carefully. More specifically, I have written insightful comments at every single line and I have also inserted arrows like this 

" <----------------- " 

which depicts the variables that I encourage the user to play with.

September 2021

Filippos Tzimkas-Dakis, Physics Dept. Univ. of Crete, Greece 
 
