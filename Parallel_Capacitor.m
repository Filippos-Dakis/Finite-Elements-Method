%% Filippos Tzimkas-Dakis  MSc. Student, UoC Physics Dept. September 2021
%
% I wrote this short because I wanted to learn some things on the FEM
% computational method. This program makes use of FEM in order to solve
% the Electrostatic problem of a Parallel Plate Capacitor. The user can
% solve the system either with Direct solver or Iterative Solvers (BicG, GMRES)
% and reveal which of them is faster. Furthermore, this program calculates
% the energy stored in the capacitor and also the capacitance with two
% different ways, FEMethod and by interpolating and integrating the result 
% obtained by solving the FEM problem.
% I encourage the user to read the program and I promise that they will
% understand every single line of it!
% 
% Do not hesitate to contact me at        
%           dakisfilippos@gmail.com    or   f.tzimkas@physics.uoc.gr


CD = cd;
addpath([CD,'\functions'])
addpath([CD,'\functions\matplotlib']);

set(groot, 'defaultAxesFontName','Times New Roman')  % changes the default FontName of Figures
%% Geometry and General properties
close all
clear all
clc
mm  = 1e-3;                       % mm multiplier
cm  = 1e-2;                       % cm multiplier

units = 1e-2;                     % change the units of the graphs (1,1e-1,1e-2,1e-3,1e-6,1e-9) !! <---------------------

w     = 3*cm;                     % (cm) plate width                  <---------------------
tt    = 2*mm;                     % (mm) plates thickness             <---------------------
d     = 1*cm;                     % (cm) distance between two plates  <---------------------
V     = 1;                        % (V) Voltage of each plate         <---------------------
V_top =  V/2;                     % (V) Upper plate Voltage           <---------------------
V_bot = -V/2;                     % (V) Lower plate Voltage           <---------------------
e0    = 8.8541878128*1e-12;       % (F/m) Vacuum permitivity
e_r   = 2.2;                      % (clean) relative permitivity      <---------------------
mu0   = 4*pi*1e-7;                % (H/m) vacuum permiability
A     = 5*w;                      % (m) X-dimension of the computational window
%------------------- GEOMETRY ---------------------------------------------
tot_rect = [3; 4; -A/2; -A/2; A/2; A/2; -A/2; A/2; A/2; -A/2];    % Geometry of the computational window
top_P = [3; 4; -w/2; -w/2; w/2; w/2;...                           % Geometry of the upper plate
    d/2; (d/2 + tt); (d/2 + tt); d/2;];
bot_P = [3; 4; -w/2; -w/2; w/2; w/2;...                           % Geometry of the lower plate
    -d/2; -(d/2 + tt); -(d/2 + tt); -d/2;];
middle = [3; 4; -w/2; -w/2; w/2; w/2;...                          % Geometry of the rectangle between the plates
    -d/2; d/2; d/2; -d/2;];                                       % because it is filled with dielectric material
gd = [tot_rect,top_P,bot_P,middle];                               % needed for matlab
ns = char('tot_rect','top_P','bot_P','middle');                   % names of the boundaries
ns = ns';                                                         % converting it to column vector
sf = 'tot_rect-top_P-bot_P+middle';                               % setting the formula for the geometry
[dl,bt] = decsg(gd,sf,ns);                                        % creates the geometry
% pdegplot(dl,'EdgeLabels','on','FaceLabels','on');   % plots the resulted geometry
% axis equal
% -------------------------------------------------------------------------
%% plot parameters for later | Dont touch them !
units_c = [1, 1e-1 1e-2 1e-3 1e-6 1e-9];       % setting the graph units 
units_C = {'m','dm','cm','mm','um','nm'};
index   = find(units_c == units);
un      = units_C{index};
x_axis  = ['X-axis (',un,')'];
y_axis  = ['Y-axis (',un,')'];
%% ---------- creates the triangularization ------------------------------
refin = 0;                           % choose the number of refinements you want  <------------------------------------

[p,e,t] = initmesh(dl);              % creates the triangularization
for ii = 1:refin
    [p,e,t] = refinemesh(dl,p,e,t);  % denses the triangular lattice
end
pdeplot(p,e,t);                      % plots the geometry with the triagualrization
title(['Lattice after ',num2str(refin),' refinement(s)'])
% -------------------------------------------------------------------------
%% Defining who node is what and applying Boundary Conditions
Nn      = size(p,2);                 % number of nodes
Ne      = size(t,2);                 % number of elements
Nd      = size(e,2);                 % number of edges - ακμές
node_id = ones(Nn,1);                % our flag in order to check if the node
%                                      is on the edge or not
X0      = zeros(Nn,1);               % initialization of the result vector

for id = 1:Nd
    if ( e(6,id )== 0 || e(7,id) == 0 )
        node_id(e(1,id)) = 0;                     % marks the specific nodes as KNOWN because they are on the edges
        node_id(e(2,id)) = 0;                     % end we have to apply the boudnary conditions
        x1 = p(1,e(1,id));                        % obtains x-coordinate for the appropriate checks
        y1 = p(2,e(1,id));                        % obtains y-coordinate for the appropriate checks
        
        if ( y1 > 0 && y1 <1.1*(d/2 + tt) && abs(x1)<= 1.1*(w/2)  )        % this IF checks if we are on the Upper Plate, if so we set the voltage
            X0(e(1,id)) = V_top;                  % Voltage = V_top.
            X0(e(2,id)) = V_top;
        elseif y1 < 0 && y1 > -1.1*(d/2 + tt) && abs(x1)<= 1.1*(w/2)       % this IF checks if we are on the Lοwer Plate, if so we set the voltage
            X0(e(1,id)) = V_bot;                  % Voltage = V_bot.
            X0(e(2,id)) = V_bot;
        else
            node_id(e(1,id)) = 1;                 % Our boundary conditions on the computational Window are
            node_id(e(2,id)) = 1;                 % von Neuman Conditions, hence we set this nodes as UNKNOWN !!
        end
    end
end
% -------------------------------------------------------------------------
counter = 0;                        % This part locates the nodes that are attributed to
for ii = 1:Nn                       % the degrees of freedom (unknown nodes) and enumerates them
    if node_id(ii) == 1             % for latter usage in the kernel;
        counter     = counter + 1;
        node_id(ii) = counter;
    end
end
% -------------------------------------------------------------------------
%% Main - Kernel (creating Stiffness and Load Matrices)
Nf  = nnz(node_id);                 % number of uknown variables
S   = spalloc(Nn,Nn,7*Nn);          % initialization of the GLOBAL Stiffness Matrix (for energy calculation)
Sff = spalloc(Nf,Nf,7*Nf);          % initialization of the large Stiff Matrix for system Solutions Sff*X = B
B   = zeros(Nf,1);                  % initialization of Load Matrix
Se  = zeros(3,3);                   % initialization of the Local Stiff Matrix
pp  = 0;
for ie = 1:Ne                       % Scans all the Elements (triangles)
    n(1:3) = t(1:3,ie);             % n(1:3) are the nodes of every elemnet
    rg     = t(4,ie);               % rg is the region of the element (subdomain number)
    x(1:3) = p(1,n(1:3));           % stores the x-coordinate
    y(1:3) = p(2,n(1:3));           % stores the y-coordinate
    
    if ( max(abs(x)) <= w/2 && max(abs(y)) <= d/2 )
        dec = e_r*e0;               % dielectric of the capacitor
    else
        dec = e0;                   % dialectric outside the capacitor
    end
    
    De   = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);  % determinant
    Ae   = abs(De/2);                                   % element area
    
    b = zeros(3,1);
    c = zeros(3,1);
    for k = 0:2
        xk = circshift(x,-k);         % Coefficients for linear
        yk = circshift(y,-k);         % triangular finite elements. For more
        b(k+1) = (yk(2)-yk(3))/De;    % information on these equations read Tsimboukis'
        c(k+1) = (xk(3)-xk(2))/De;    % complementary books "Notes on Computational Electromagnetics" pages 83-95
    end                               % and "Energy Methods for E/M Fields" pages 287-308
    for i = 1:3                       % calculates the 3X3 Stiff Local Matrix
        for j = 1:3                   %
            Se(i,j)      = ( b(i)*b(j) + c(i)*c(j) )*dec*Ae;  % Local Stiff Matrix
            S(n(i),n(j)) = S(n(i),n(j)) + Se(i,j);            % GLOBAL STIFFNESS MATRIX for ENERGY Calculation
            if ( node_id(n(i)) ~= 0 )                         % node n(i) is Unknown
                if ( node_id(n(j)) ~= 0 )                     % node n(j) is Unknown
                    Sff(node_id(n(i)),node_id(n(j))) = Sff(node_id(n(i)),node_id(n(j))) + Se(i,j);
                else
                    % calcultes the Load Matrix
                    B(node_id(n(i))) =  - Se(i,j)*X0(n(j)) + B(node_id(n(i)));
                end
            end
        end
    end
end
%--------------------------------------------------------------------------
%% ----- Direct Solver ----------------------------------------------------
tic
X = Sff\B;                       % solves the system directly  Ax = B
toc
% -------------------------------------------------------------------------
%% -------- iterative solver: "Biconjugate gradient" -> Ax = B ------------
% [x,flag] = bicg(A,b,tol,maxit)
% If tol is [], then bicg uses the default, tol=1e-6
% maxit ----> maximum number of iterations
% flag = 0 -> bicg converged to the desired tolerance tol within maxit iterations
% flag = 1 -> bicg iterated maxit times but did not converge.
bicg_ON_OFF = 0;                  % set 1 if you want to solve the system with    <---------------------
if bicg_ON_OFF                    % BICG iterative solver
    flag1 = 1;
    ii    = 0;
    while flag1 == 1              % We use this While-loop to find out the least number of iterations needed
        ii         = ii + 1;      % so as Iterative Solver  bicg()  to converge;
        [~, flag1] = bicg(Sff,B,[],ii);
    end
    if ~flag1
        tic
        [X1, flag1] = bicg(Sff,B,[],ii);
        toc
        fprintf('\n')
        cprintf('*cyan',' IT solver BicG converged at %d iterations. \n',ii);
    end
end
% -------------------------------------------------------------------------
%% --------- iterative solver: "GMRES" -> Ax = b --------------------------
% [x,flag] = gmres(A,b,restart,tol)
% If tol is [], then bicg uses the default, tol=1e-6
% maxit ----> maximum number of iterations bounded by SIZE(A,1)
% flag = 0 -> bicg converged to the desired tolerance tol within maxit iterations
% flag = 1 -> bicg iterated maxit times but did not converge.
gmres_ON_OFF = 0;                % set 1 if you want to solve the system with    <---------------------
if gmres_ON_OFF                  % GMRES iterative solver
    flag2 = 1;
    ii    = 0;
    while flag2 == 1             % We use this While-loop to find out the least number of iterations needed
        ii         = ii + 1;     % so as Iterative Solver   gmres()   to converge;
        [~, flag2] = gmres(Sff,B,[],[],ii);
    end
    if ~flag2
        tic
        [X2, flag2] = gmres(Sff,B,[],[],ii);
        toc
        fprintf('\n')
        cprintf('*yellow',' IT solver GMRES converged at %d iterations. \n',ii);
    end
end
% -------------------------------------------------------------------------
%% Forming the total vector
for ii = 1:length(X0)            % completing the result vector after the system solution
    if node_id(ii) > 0           % by creatinga total vector including the boundaries
        X0(ii) = X(node_id(ii)); % choose X, X1 or X2
    end
end
%% Converts the geometry into "units"meter   (m,cm,mm etc)
gd = [tot_rect,top_P,bot_P];                             % needed for matlab
ns = char('tot_rect','top_P','bot_P');                   % names of the boundaries
ns = ns';                                                % converting it to column vector
sf = 'tot_rect-top_P-bot_P';                             % setting the formula for the geometry
[dl2,bt] = decsg(gd,sf,ns);                              % for better interpretation
dl_plot = dl2;
dl_plot(2:5,1:end) = dl2(2:5,1:end)/units;
%% ------------------  PLOTS  ---------------------------------------------
figure
subplot(1,2,2)
[Ex,Ey] = pdegrad(p,t,-X0);                        % calculates the Electric field
pltE = pdeplot(p/units,e,t,'FlowData',[Ex;Ey]);    % plots the Electric Field
pltE.LineWidth = 1;
pltE.Color =  [1 1 0];
hold on,box on
set(gca,'Color','k')                               % sets the background color black
h = pdegplot(dl_plot);                             % plots the resulted geometry
h.Color = [1 1 1];                                 % with white edges
box on, axis equal
xlabel(x_axis),ylabel(y_axis)
title('Electric Field'),hold off;


subplot(1,2,1)
pdeplot(p/units,e,t,'XYData',X0,'Contour','on');     % plots the Potential in colormap
view(2)
set(gca,'Color','k')
pos = [0.1    0.295     0.3329    0.4438];           % sets a specific positio for the plot
set(gca,'position',pos)
hold on, box on
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.Color = [1 1 1];                                   % with white color
colormap(viridis)                                    % hot,parula,cool,summer,gray, winter,
xlabel(x_axis), ylabel(y_axis)
axis equal
title('Electric Potential'),hold off;
%%   Energy and Capacitance via FEM
% Calculates the energy inside the capacitor by taking into account only
% the elements between the two plates
W_e = 0;
for ie = 1:Ne                       % Scans all the Elements (triangles)
    n(1:3) = t(1:3,ie);             % n(1:3) are the nodes of every elemnents
    rg     = t(4,ie);               % rg is the region of the element (subdomain number)
    x(1:3) = p(1,n(1:3));           % stores the x-coordinate
    y(1:3) = p(2,n(1:3));           % stores the y-coordinate
    if ( max(abs(x(1:3)))<= w/2 && max(abs(y(1:3))) <= d/2 )    % Considering only the elements with coordinates
        %                                                         between the two plates |x|<=w/2 and |y|<=d/2
        dec = e_r*e0;                                           % dielectric of the capacitor
        
        De   = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);      % determinant
        Ae   = abs(De/2);                                       % element area
        
        b = zeros(3,1);
        c = zeros(3,1);
        for k = 0:2
            xk = circshift(x,-k);         % Coefficients for linear
            yk = circshift(y,-k);         % triangular finite elements. For more
            b(k+1) = (yk(2)-yk(3))/De;    % information on these equations read Tsimboukis'
            c(k+1) = (xk(3)-xk(2))/De;    % complementary books "Notes on Computational Electromagnetics" pages 83-95
        end                               % and "Energy Methods for E/M Fields" pages 287-308
        for i = 1:3                    % calculates the 3X3 Stiff Local Matrix
            for j = 1:3                %
                Se(i,j) = ( b(i)*b(j) + c(i)*c(j) )*dec*Ae;  % Local Stiff Matrix
                W_e     = W_e + 1/2 * X0(n(i)) * Se(i,j) * X0(n(j));
                
            end
        end
    end
end

C_fem = 2*W_e/V^2;                            % calculates capacitance of between the two plates
fprintf('\n')                             
cprintf('*white',' C_FEM = %.3f pf \n',C_fem*1e12);   % informs the user for the result
%% ---- Data Interpolation ------------------------------------------------
% Interpolates the Electric potential and then calculates the electric field.
F       = pdeInterpolant(p,t,X0);                     % preparation for interpolation command
XX1     = linspace(-A/2,A/2,1000);                    % points along X-axis
YY1     = linspace(-A/2,A/2,1000);                    % points along Y-axis
dX1     = abs(XX1(1) - XX1(2));                       % x-step
dY1     = abs(YY1(1) - YY1(2));                       % y-step
[X1,Y1] = meshgrid(XX1,YY1);                          % grid generation
uOut    = evaluate(F,X1,Y1);                          % interpolation command execution
Z       = reshape(uOut,size(X1,1),size(Y1,1));        % from vector to array
[Fx,Fy] = gradient(Z);                                % calcultes the Electric Field
Fx      = -Fx/dX1;                                    % inserting the x-step to derivative
Fy      = -Fy/dY1;                                    % inserting the y-step to derivative
Z2      = sqrt( abs(Fx).^2 + abs(Fy).^2 );            % calculates the field amplitude

figure
surf(X1/units,Y1/units,Z2,'EdgeColor','none');        % plots field amplitude
view(2)                                               % 2D vies
colormap(viridis),colorbar
hcb = colorbar;
hcb.Label.String = 'E ( V/m )';
axis equal, hold on
set(gca,'Color','k')                                 % sets the background Black
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                   % with white color
xlabel(x_axis),ylabel(y_axis)
title('Electric Field in Capacitor')

%% Energy and Capacitance through Intepolation (post processing) ----------
W_int = Energy_Interp(p,t,X0,e_r*e0,w/2,d/2,500);            % Energy between the capacitor plates with interpolation
% W_int = Energy_Interp(p,t,X0,epsilon,x_max,y_max,points)
C_int    = 2*W_int/(V_top-V_bot)^2;                           % capacitance through interpolation results
C_theor  = e0*e_r*w/d;                                        % simplified theoretical result
C_int/C_theor;                                                % comparing the 2 capacitances
fprintf('\n')
cprintf('*white',' C_int = %.3f pf \n',C_int*1e12);           % informs the user for the result
fprintf('\n')
cprintf('*white',' C_theor = %.3f pf \n',C_theor*1e12);       % informs the user for the result
% -------------------------------------------------------------------------