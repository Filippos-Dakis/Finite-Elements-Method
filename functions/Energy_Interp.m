function W_interp = Energy_Interp(p,t,X0,epsilon,x_max,y_max,points)
%% Filippos Tzimkas-Dakis  MSc. Student, UoC Physics Dept. September 2021
%
% This function calculates the energy for an Electrostatic proflem.
% p,t,X0 is the grid parameters and X0 the electric potential φ calculated
%        with FEMethod
% epsilon = e_r*e0    (F/m) dielectric constant 
% x_max,y_max the --> [-x_max,x_max]&[-y_max,y_max] dimensions in which we 
%                     want to interpolate φ potential and then calculate 
%                     the El. Field
% points = number of interpolation points in each direction
%
% Do not hesitate to conatct me at    
%           dakisfilippos@gmail.com    or   f.tzimkas@physics.uoc.gr
%

x_max   = abs(x_max);
y_max   = abs(y_max);

F       = pdeInterpolant(p,t,X0);                     % preparation for interpolation command
XX1     = linspace(-x_max,x_max,points);          % points along X-axis
YY1     = linspace(-y_max,y_max,points);          % points along Y-axis
dX1     = abs(XX1(1) - XX1(2));                       % x-step
dY1     = abs(YY1(1) - YY1(2));                       % y-step
[X1,Y1] = meshgrid(XX1,YY1);                          % grid generation
uOut    = evaluate(F,X1,Y1);                          % interpolation command execution
Z       = reshape(uOut,size(X1,1),size(Y1,1));        % from vector to array
[Fx,Fy] = gradient(Z);                                % calcultes the Electric Field
Fx      = -Fx/dX1;                                    % inserting the x-step to derivative
Fy      = -Fy/dY1;                                    % inserting the y-step to derivative
Z2      = sqrt( abs(Fx).^2 + abs(Fy).^2 );            % calculates the field amplitude

% Energy = (1/2)*epsilon* \int\int |E|^2 dS 
W_interp = (1/2)*epsilon*sum( abs(Z2(~isnan(Z2))).^2)*dX1*dY1;  % energy inside the capacitor
end

