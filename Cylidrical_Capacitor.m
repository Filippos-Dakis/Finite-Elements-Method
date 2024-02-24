%% Filippos Tzimkas-Dakis  MSc. Student, UoC Physics Dept. September 2021
%
% I wrote this short because I wanted to learn some things on the FEM
% computational method. There is a cylindrical capacitor (infinite length),
% outer radius R_out and inner radius R_in. The outer cylinder has potential
% V = V_out volt while the the inner cylinder has electric potential V = V_in .
% For educational purpose I solve the problem with three(3) different
% solvers, one direct solver (D.S.) and two iterative solvers (It.S.) .
% There are plenty of comments so I strongly encourage the user to read the
% code and I am confident that they will understand every single step. 
% 
% Do not hesitate to contact me at        
%           dakisfilippos@gmail.com    or   f.tzimkas@physics.uoc.gr


CD = cd;
addpath([CD,'\functions'])
addpath([CD,'\functions\matplotlib']);

set(groot, 'defaultAxesFontName','Times New Roman')  % changes the default FontName of Figures
%% 1 General & Geometry
close all
clear all
clc
mm = 1e-3;                          % set the length units of the radii you will use 
cm = 1e-2;                          %  m->1, cm->1e-2, mm-> 1e-3 etc

units = 1e-3;                       % change the units of the graphs (1,1e-1,1e-2,1e-3,1e-6,1e-9) !!  <---------------------

R_out = 3.5/2;                      % (m_units) outer radius            <---------------------
R_in  = 1/2;                        % (m_units) inner radius            <---------------------
R_out = R_out * mm;                 % setting the dimensions correctly  <---------------------
R_in  = R_in  * mm;                 % setting the dimensions correctly  <---------------------
e0    = 8.8541878128*1e-12;         % (F/m) Vacuum permitivity 
e_r   = 1;                          % (clean) relative permitivity
mu0   = 4*pi*1e-7;                  % (H/m) vacuum permiability
Z0    = 50;                         % (Ohm) Characteristic Impedance of the "cable"  <---------------------
coaxial_cable = 1;                  % set it 1 if you want the cable to have         <---------------------
if coaxial_cable                    % Char.Imp. Z = Z0 Ohm  (usualy Z0 = 50).
    R_in = R_out * exp(-sqrt(e_r)*2*pi*Z0/sqrt(mu0/e0));    % (m) inner radius
    % The above relation comes from the known, for coaxial cables, relation
    % Z0 = sqrt(L/C) = 1/(2*pi) * sqrt(mu0/e0) * sqrt(1/e_r) * ln(R_out/R_in)    
end

V_out = 0;             % Electirc potential on the outer cylinder   <---------------------
V_in  = 1;             % Electirc potential on the inner cylinder   <---------------------

% creating the geometry. search for 'decsg' in MATLAB documentation
C1 = [1 0 0 R_out]';   % outter circle; column vector
C2 = [1 0 0 R_in]';    % inner circle;  column vector
gd = [C1,C2];          % needed for matlab
ns = char('C1','C2');  % names of the boundaries
ns = ns';              % converting it to column vector
sf = 'C1-C2';          % setting the formula for the geometry

[dl,bt] = decsg(gd,sf,ns);  % creates the geometry

% pdegplot(dl,'EdgeLabels','on','FaceLabels','on');   % plots the resulted geometry
% xlim([-R_out,R_out])
% axis equal
%% plot parameters for later | Dont touch them !
units_c = [1, 1e-1 1e-2 1e-3 1e-6 1e-9];       % setting the graph units 
units_C = {'m','dm','cm','mm','um','nm'};
index   = find(units_c == units);
un      = units_C{index};
x_axis  = ['X-axis (',un,')'];
y_axis  = ['Y-axis (',un,')'];
%% -------------- creating the mesh  --------------------------------------
refin   = 2;                         % choose how many refinments you want  <---------------------
[p,e,t] = initmesh(dl);              % creates the triangularization
for ii = 1:refin
    [p,e,t] = refinemesh(dl,p,e,t);  % denses the triangular lattice
end
pdeplot(p/units,e,t);                      % plots the geometry with the triagualrization
xlabel(x_axis)
ylabel(y_axis)
title(['Lattice after ',num2str(refin),' refinement(s)'])

%% 2 "Re-number" the unknown nodes
Nn = size(p,2);                                   % number of nodes
Ne = size(t,2);                                   % number of elements
Nd = size(e,2);                                   % number of edges - ακμές
node_id = ones(Nn,1);                             % our flag in order to check if the node
%                                                    is on the edge or not
X0      = zeros(Nn,1);                            % initialization of the result vector
                                                  % In this part we apply the boundary conditions.
for id = 1:Nd
    if ( e(6,id )== 0 || e(7,id) == 0 )
        node_id(e(1,id)) = 0;                     % marks the specific nodes as KNOWN because they are on the edges
        node_id(e(2,id)) = 0;                     % end we have to apply the boudnary conditions 
        radd = sqrt( p(1,e(1,id))^2 +  p(2,e(1,id))^2 );   % calculates the radius of the node  R = sqrt( x^2 + y^2 );
        if radd > ( R_in + (R_out - R_in)/2 )     % this if checks if we are at the edge of the outter or the inner
            X0(e(1,id)) = V_out;                  % cylinder. If we are at the outter edge we set the Voltage =  V_out
            X0(e(2,id)) = V_out;                  % All the others remain zero, both the unknow and the Inner edge.
        else
            X0(e(1,id)) = V_in;                   % sets the Voltage on the inner cylinder Voltage = V_in
            X0(e(2,id)) = V_in;
        end
    end
end
% -------------------------------------------------------------------------
counter = 0;                       % This part locates the nodes that are attributed to 
for ii = 1:Nn                      % the degrees of freedom (unknown nodes) and enumerates them
    if node_id(ii) == 1            % for latter usage in the kernel;
        counter = counter + 1;
        node_id(ii) = counter;
    end
end
% -------------------------------------------------------------------------
%% 3 Kernel
Nf  = nnz(node_id);                % number of uknown variables
S   = spalloc(Nn,Nn,7*Nn);         % initialization of the GLOBAL Stiffness Matrix (for energy calculation) 
Sff = spalloc(Nf,Nf,7*Nf);         % initialization of the large Stiff Matrix for system Solutions Sff*X = B
B   = zeros(Nf,1);                 % initialization of Load Matrix
Se  = zeros(3,3);                  % initialization of the Local Stiff Matrix
pp  = 0;
for ie = 1:Ne                      % Scans all the Elements (triangles)
    n(1:3) = t(1:3,ie);            % n(1:3) are the nodes of every elemnet
    rg     = t(4,ie);              % rg is the region of the element (subdomain number)
    x(1:3) = p(1,n(1:3));          % stores the x-coordinate
    y(1:3) = p(2,n(1:3));          % stores the y-coordinate
    De   = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);  % determinant 
    Ae   = abs(De/2);              % element area

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
            Se(i,j) = ( b(i)*b(j) + c(i)*c(j) )*e_r*Ae;  % Local Stiff Matrix
            S(n(i),n(j)) = S(n(i),n(j)) + Se(i,j);       % GLOBAL STIFFNESS MATRIX for ENERGY Calculation
            if ( node_id(n(i)) ~= 0 )                    % node n(i) is Unknown
                if ( node_id(n(j)) ~= 0 )                % node n(j) is Unknown
                    Sff(node_id(n(i)),node_id(n(j))) = Sff(node_id(n(i)),node_id(n(j))) + Se(i,j);
                else
                    % calcultes the Load Matrix 
                    B(node_id(n(i))) =  - Se(i,j)*X0(n(j)) + B(node_id(n(i)));  
                end
            end
        end
    end
end
% ----- Direct Solver --------------
tic
X = Sff\B;                       % solves the system directly  Ax = B
toc
% ----------------------------------
% ---------iterative solver: "Biconjugate gradient" -> Ax = B -------------
% [x,flag] = bicg(A,b,tol,maxit)
% If tol is [], then bicg uses the default, tol=1e-6
% maxit ----> maximum number of iterations
% flag = 0 -> bicg converged to the desired tolerance tol within maxit iterations
% flag = 1 -> bicg iterated maxit times but did not converge.
flag1 = 1;                          %                             <---------------------
ii    = 0;
while flag1 == 1                    % We use this While-loop to find out the least number of iterations needed
    ii         = ii + 1;            % so as Iterative Solver  bicg()  to converge;
    [~, flag1] = bicg(Sff,B,[],ii); 
end
if ~flag1
    tic
    [X1, flag1] = bicg(Sff,B,[],ii);
    toc
    fprintf('\n')
    cprintf('*cyan',' IT solver BicG converged at %d iterations. \n',ii);
end
% -------------------------------------------------------------------------
% --------- iterative solver: "GMRES" -> Ax = b ---------------------------
% [x,flag] = gmres(A,b,restart,tol)
% 
% If tol is [], then bicg uses the default, tol=1e-6
% maxit ----> maximum number of iterations bounded by SIZE(A,1)
% flag = 0 -> bicg converged to the desired tolerance tol within maxit iterations
% flag = 1 -> bicg iterated maxit times but did not converge.
flag2 = 1;                              %                       <---------------------
ii    = 0;
while flag2 == 1                        % We use this While-loop to find out the least number of iterations needed
    ii         = ii + 1;                % so as Iterative Solver   gmres()   to converge;
    [~, flag2] = gmres(Sff,B,[],[],ii); 
end
if ~flag2
    tic
    [X2, flag2] = gmres(Sff,B,[],[],ii);
    toc
    fprintf('\n')
    cprintf('*yellow',' IT solver GMRES converged at %d iterations. \n',ii);
end
% -------------------------------------------------------------------------

for ii = 1:length(X0)            % completing the result vector after the system solution
    if node_id(ii) > 0           % by creatinga total vector including the boundaries
        X0(ii) = X(node_id(ii)); % choose X, X1 or X2 
    end
end

%%   Energy and Capacitance 
W      = (pi*(V_in/log(R_out/R_in)) ^2)*log(R_out/R_in);   % Theoretical Energy
W_FEM  = (1/2) * (X0')*S*X0;                               % Energy calculated with FEM 
WW     = W_FEM/W ;                                         % Comparing the two energies
C_FEM  = 2*W_FEM/(V_in^2);                                 % Capacitance in units of \epsilon_0 (F/m) calculated with FEM results
C      = 2*pi*e_r/log(R_out/R_in);                         % Capacitance in units of \epsilon_0 (F/m) calculated with analytical expression
CC     = C_FEM/C;                                          % comparing the two capacitances
%% Converts the geometry into "units"meter   (m,cm,mm etc)
dl_plot = dl;
dl_plot(2:end,1:end) = dl(2:end,1:end)/units;
%% ------------------  PLOTS  --------------------------------------------- 
x_max = max(p(1,:));                               % max x-position
x_min = min(p(1,:));                               % min x-position
y_max = max(p(2,:));                               % max y-position
y_min = min(p(2,:));                               % min y-position

figure
subplot(1,2,2)
[Ex,Ey] = pdegrad(p,t,-X0);                        % calculates the Electric field
pltE = pdeplot(p/units,e,t,'FlowData',[Ex;Ey]);    % plots the Electric Field
pltE.LineWidth = 1;
set(gca,'Color','k')                               % sets the background color black
hold on                                             
h = pdegplot(dl_plot);                             % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                 % with white edges
box on, axis equal
xlim([x_min x_max]*1.1/units)                      % xlimits
ylim([y_min y_max]*1.1/units)                      % ylimits
xlabel(x_axis)
ylabel(y_axis)
title('Electric Field')
hold off


subplot(1,2,1)
pdeplot(p/units,e,t,'XYData',X0);                  % plots the Potential in a colormap
set(gca,'Color','k')
pos = [0.1    0.295     0.3329    0.4438];         % sets a specific positio for the plot 
set(gca,'position',pos)
hold on
h = pdegplot(dl_plot);                             % plots the resulted geometry
h.Color = [1 1 1];                                 % with white color
box on
colormap(viridis)                                  % hot,parula,cool,summer,gray, winter,
xlabel(x_axis)
ylabel(y_axis)
axis equal
xlim([x_min x_max]*1.1/units)                      % xlimits
ylim([y_min y_max]*1.1/units)                      % ylimits
title('Electric Potential')
hold off
sg = sgtitle({'Elec. Potential $\Phi$ \& Elec. Field $\rm {\vec{E}}$'},'Interpreter','latex');
sg.FontName    = 'Times';
%% ----------------- Data Interpolation for Energy Calculation ------------
% 
M       = 500;                          % number of points for interpolation   <---------------------
N       = M;                            % number of points for interpolation
R1      = R_in;                         % inner radius 
R2      = R_out;                        % outer radius
nR      = linspace(R1,R2,M) ;           % setting the nodes along rho coordinate
nT      = linspace(0,2*pi,N) ;          % setting the nodes along phi coordinate
                                        % nT = pi/180*(0:NT:theta) ;
dphi    = abs(nT(1) - nT(2));
dR      = abs(nR(1) - nR(2));
[R, T]  = meshgrid(nR,nT) ;             % creates the mesh
X1      = R.*cos(T);                    % Convert grid to cartesian coordintes
Y1      = R.*sin(T);                    % Convert grid to cartesian coordintes
% ---------- Data interpolation to cartesian coordinates ------------------
F       = pdeInterpolant(p,t,X0);
uOut    = evaluate(F,X1,Y1);
Z       = reshape(uOut,size(X1,1),size(Y1,1));
% -------------------------------------------------------------------------
% ----------- plotting interpolated data ----------------------------------
figure
subplot(1,2,1)
surf(X1/units,Y1/units,Z,'EdgeColor','none')    % Electric Potential 
view(2)
colormap(viridis)
hcb = colorbar;
hcb.Label.String = '\Phi ( V )';
axis equal
set(gca,'Color','k')
hold on
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                   % with white color
xlabel(x_axis)
ylabel(y_axis)
xlim([x_min x_max]*1.1/units)                        % xlimits
ylim([y_min y_max]*1.1/units)                        % ylimits
title('Potential')

%        ----------------------------------------
[dFr,dFt] = gradient(Z,nR,nT);                       % calculates the Electromagnetic Field
dFt       = -dFt./(repmat(nR(:)',length(nR),1));
dFr       = -dFr;
Z2        = sqrt( abs(dFr).^2 + abs(dFt).^2 );       % calculates the Electric Field Amplitude
subplot(1,2,2)
surf(X1/units,Y1/units,Z2,'EdgeColor','none')
view(2)
colormap(viridis)
hcb = colorbar;
hcb.Label.String = 'E ( V/m )';
axis equal, hold on
set(gca,'Color','k')
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                   % with white color
xlabel(x_axis)
ylabel(y_axis)
xlim([x_min x_max]*1.1/units)                        % xlimits
ylim([y_min y_max]*1.1/units)                        % ylimits
title('Electric Field')

sg = sgtitle({'Elec. Potential $\Phi$ \& Elec. Field $\rm {|\vec{E}|}$'},'Interpreter','latex');
sg.FontName = 'Times';
%% ------------ numerical Polar integration   -----------------------------
ZZ = 0;                                        
for ii = 1:size(Z2,1)                                    % ZZ2 Rows correspond to constant angle \phi and the radius spanning
    cc = sum(~isnan(Z2(ii,:)));                          % from R_in to R_out. Hence, |E|^2 * ρ = (Z2(ii,1:cc).^2).*nR(1:cc)
    ZZ = ZZ + sum( (Z2(ii,1:cc).^2).*nR(1:cc) )*dR*dphi; % We also have to multiply with dR and dPhi because of the integration
end

W_interp = (1/2)*ZZ;                                     % Energy calculated via the numerical polar integration
C_interp = 2*W_interp/((V_in- V_out)^2) ;                % The corresponding Capacitance

%% phi -> theoretical approach
K  = V_in/log(R_in/R_out);
pp = linspace(R_in, R_out, 200);
v  = K*log(pp./R_out);

figure 
plot(pp/units,v,'Color','b','LineWidth',1.5)
xlim([R_in,R_out]/units)
xlabel(['\rho (',un,')'])
ylabel('\phi (V)')
title('$ \phi(\rho) = -V\frac{ln(\rho/b)}{ln(b/a)} $','Interpreter','latex')

line([pp(100),pp(100)]/units,[0,v(100)],'Color','red','LineStyle','--')
line([0,pp(100)]/units,[v(100),v(100)],'Color','red','LineStyle','--')
text(pp(100)/units,v(100),'$ \leftarrow \phi=0.4,\,  \rho=(a+b)/2 $','Interpreter','latex')








