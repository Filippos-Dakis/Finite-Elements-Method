%% Filippos Tzimkas-Dakis MSc. Student, UoC Physics Dept. September 2021
%
% I wrote this short because I wanted to learn some things on the FEM
% computational method. In this program I am solving the scattering
% problem of an infinite(along z-axis) PEC cylinder. I assume an incident
% wave with TM polarization (i.e. E||z-axis) which propagates parallel to
% x-axis (from negative x to positive x). In order to make solution
% convenient I use PML layers around the computational window. We test the
% consistency of the result by also ploting the analytical solution given
% by an infite expansion of second order Hankel fucntions. Also, in order
% to understand the geometry of the problem, have a look at the file named
% "Scattering_Problem.pdf".
% There are plenty of comments so I strongly encourage the user to read the
% code and I am confident that they will understand every single step. 
% 
% Do not hesitate to contact me at        
%           dakisfilippos@gmail.com    or   f.tzimkas@physics.uoc.gr

CD = cd;
addpath([CD,'\functions'])
addpath([CD,'\functions\matplotlib']);

set(groot, 'defaultAxesFontName','Times New Roman')  % changes the default FontName of Figures
%%
close all, clear all,clc
xC     = 0;                          % x-axis cylinder center point
yC     = 0;                          % y-axis cylinder center point
e_r    = 1;                          % (clean) relative permitivity     <---------------------
e0     = 8.8541878128*1e-12;         % (F/m) Vacuum permitivity;
mu0    = 4*pi*1e-7;                  % (H/m) vacuum permiability
freq   = 3e8;                        % (Hz)incident field frequency     <---------------------
c0     = 3e8;                        % (m/s) speed of light   1/sqrt(e_r*e0*mu0);
lambda = c0/freq;                    % (1/m) wavelength in the computational window
w      = 1*lambda;                   % (m) distance between scatterer and PML   <---------------------
Eo     = 1;                          % (V/m) incident field amplitude           <---------------------

Cyl_Rad= lambda/4 ;                  % (m) Cylinder's radius                    <---------------------

d      = lambda/4;                   % (m) width of PML layer
omega  = 2*pi*freq;                  % (rad/s) angular frequency
k0     = omega*sqrt(e_r*e0*mu0);     % (rad/m) wavenumber
pml_ref= 1e-6;                       % PML reflectivity for normal incident wave   <---------------------
beta   = -log(pml_ref)/(2*k0*d);     % \beta PML coefficient
alpha  = 1-1i*beta;                  % \alph PML coefficient
llx    = [1/alpha alpha alpha];      % pml vector along x-direction
lly    = [alpha  1/alpha alpha];     % pml vector along y-direction
llxy   = [1 1 alpha^2];              % pml at the the corners of the window

%% creating the differenct regions
x = zeros(9,4);
y = zeros(9,4);

for ii = 1:2                                % creates the computational window
    x(1,ii)     =  ((-1)^ii)*(w + Cyl_Rad);
    y(1,ii)     =  (w + Cyl_Rad);
    x(1,ii + 2) =  ((-1)^ii)*(w + Cyl_Rad);
    y(1,ii + 2) = -(w + Cyl_Rad);
end

x_max = max(x(1,:));                 % max x-position (computational window)
x_min = min(x(1,:));                 % min x-position (computational window)
y_max = max(y(1,:));                 % max y-position (computational window)
y_min = min(y(1,:));                 % min y-position (computational window)

for jj = 2:3                                % creates the Left and Right PML windows
    for ii = 1:2
        x(jj,ii)    = (-1)^(jj - 1) *(w + Cyl_Rad + d);
        x(jj,ii +2) = (-1)^(jj - 1) *(w + Cyl_Rad);
        y(jj,ii)    = (-1)^(ii - 1) *(w + Cyl_Rad);
        y(jj,ii +2) = ((-1)^ii)     *(w + Cyl_Rad);
    end
end
for jj = 4:5                                % creates the Upper and Lower PML windows
    for ii = 1:2
        y(jj,ii)    = (-1)^(jj - 1) *(w + Cyl_Rad + d);
        y(jj,ii +2) = (-1)^(jj - 1) *(w + Cyl_Rad);
        x(jj,ii)    = ((-1)^ii) *(w + Cyl_Rad);
        x(jj,ii +2) = ((-1)^ii) *(w + Cyl_Rad);
    end
end
for jj =6:7                                  % creates the Lower XY overlapping PML regions
    for ii =1:2
        x(jj,ii)    = ((-1)^(jj-1)) * (w + Cyl_Rad + d);
        x(jj,ii +2) = ((-1)^(jj-1)) * (w + Cyl_Rad);
        y(jj,ii)    = -(w + Cyl_Rad + d*(ii-1));
        y(jj,ii +2) = -(w + Cyl_Rad + d*abs(ii-2));
    end
end
for jj =8:9                                  % creates the Upper XY overlapping PML regions
    for ii =1:2
        x(jj,ii)    = ((-1)^jj) * (w + Cyl_Rad + d);
        x(jj,ii +2) = ((-1)^jj) * (w + Cyl_Rad);
        y(jj,ii)    = (w + Cyl_Rad + d*(ii-1));
        y(jj,ii +2) = (w + Cyl_Rad + d*abs(ii-2));
    end
end
c1 = 3*ones(1,9);
c2 = 4*ones(1,9);
gd=[c1         1;                     % Creates the matrix for the geometry
    c2        xC;                     % Last column consists of the PEC cylinder
    x(:,1)'   yC;
    x(:,2)' Cyl_Rad;
    x(:,4)'    1;
    x(:,3)'    1;
    y(:,1)'    1;
    y(:,2)'    1;
    y(:,3)'    1;
    y(:,4)'    1];
sf = 'R1-R0+R2+R3+R4+R5+R6+R7+R8+R9';                             % setting the formula for the geometry
ns = char('R1','R2','R3','R4','R5','R6','R7','R8','R9','R0');     % names of the boundaries
ns = ns';                                                         % converting it to column vector
[dl,bt] = decsg(gd,sf,ns);                                        % creates the geometry
figure 
subplot(1,2,1)
pdegplot(dl,'FaceLabels','on');              % plot the resulted geometry
axis equal, axis tight
xlabel('X-axis (m)'),ylabel('Y-axis (m)')    % XY-label
title('Geometry & Regions'), hold off
%% ---------- creates the triangularization ------------------------------
refin = 3;                           % choose the number of refinements you want  <---------------------
[p,e,t] = initmesh(dl);              % creates the triangularization
for ii = 1:refin
    [p,e,t] = refinemesh(dl,p,e,t);  % denses the triangular lattice
end
subplot(1,2,2)
pdeplot(p/lambda,e,t);               % plot the geometry with the triagualrization
axis equal, axis tight
xlabel('X-axis (\lambda)'),ylabel('Y-axis (\lambda)')   % XY-labels
title(['Lattice after ',num2str(refin),' refinement(s)']),hold off

%% Preparations and initializations for the MAIN part 
Nn      = size(p,2);   % total number of nodes
Ne      = size(t,2);   % total number of elements
Nd      = size(e,2);   % number of edges - ακμές
node_id = ones(Nn,1);  % to identify if a node is KNOWN or UNKNOWN
Ei      = zeros(Nn,1); % incident field
Es      = zeros(Nn,1); % Scattered Field

% Calculte the incident field Ei in the main computational window
for ie = 1:Ne
    n(1:3) = t(1:3,ie);
    rg     = t(4,ie);
    x(1:3) = p(1,n(1:3));
    y(1:3) = p(2,n(1:3));
    if rg == 9       % the largerest area is alaways at the bottom line
        for ii = 1:3
            Ei(n(ii)) = Eo * exp(-1i*k0*(x(ii)));
        end
    end
end

% Boundary conditions
for ii = 1:Nd
    if (e(6,ii) == 0 || e(7,ii) == 0 )
        for jj = 1:2
            % the following IF checks if we are on the boundary of the PEC cylinder or
            % the total window
            if ( abs( p(1,e(jj,ii)) ) < (Cyl_Rad + w )  &&  abs( p(2,e(jj,ii)) ) < (Cyl_Rad + w) )
                node_id(e(jj,ii)) = 0;                % defines the node as known
                Es(e(jj,ii))      = - Ei(e(jj,ii));   % defines the Dirichlet condition on the PEC cylider
            end
        end
    end
end
% Re-enumerating the UNKNOWNS !
counter = 0;
index   = zeros(Nn,1);
for ii = 1:Nn
    if node_id(ii) == 1
        counter = counter + 1;
        index(ii) = counter;     % index for the unknowns
    end
end
%% Main - Kernel Calculation of Stiffness and Mass Matrices ------------
Nf  = nnz(node_id);                 % number of uknown variables
Sff = spalloc(Nf,Nf,7*Nf);          % initialization of the large Stiff Matrix for system Solutions Sff*X = B
B   = zeros(Nf,1);
Aet = spalloc(Nf,Nf,7*Nf);


muxx(1) = mu0;         % computational window, the area where field is being
muyy(1) = mu0;         % calculated
ezz(1)  = e_r*e0;
%
muxx(2) = mu0*llx(1);  % left and right X regions with PML
muyy(2) = mu0*llx(2);
ezz(2)  = e0* llx(3);
%
muxx(3) = mu0*lly(1);  % top and bottom Y regions with PML
muyy(3) = mu0*lly(2);
ezz(3)  = e0* lly(3);
%
muxx(4) = mu0*llxy(1); % overlapping xy regions with PML
muyy(4) = mu0*llxy(2);
ezz(4)  = e0* llxy(3);

% -------- MAIN LOOP ------------------------------------------------------
S  = zeros(3,3);
Te = zeros(3,3);
T  = (1/12) * ones(3);                  % local Mass Matrix is standard !
T  = T - (1/12)*eye(3) + (1/6)*eye(3);  % Te = Ae* [1/6    1/12   1/12;
%                                                   1/12   1/6    1/12;
%                                                   1/12   1/12    1/6 ];
tic
for ie = 1:Ne
    n(1:3) = t(1:3,ie);            % n(1:3) are the nodes of every elemnet
    rg     = t(4,ie);              % rg is the region of the element (subdomain number)
    x(1:3) = p(1,n(1:3));          % stores the x-coordinate
    y(1:3) = p(2,n(1:3));          % stores the y-coordinate
    De     = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);  % determinant
    Ae     = abs(De/2);            % element area
    
    b(1) = ( y(2) -y(3) ) / De ;    % Coefficients for linea
    b(2) = ( y(3) -y(1) ) / De ;    % triangular finite elements. For more
    b(3) = ( y(1) -y(2) ) / De ;    % information on these equations read Tsimboukis'
    c(1) = ( x(3) -x(2) ) / De ;    % complementary books "Notes on Computational Electromagnetics" pages 83-95
    c(2) = ( x(1) -x(3) ) / De ;    % and "Energy Methods for E/M Fields" pages 287-308
    c(3) = ( x(2) -x(1) ) / De ;
    % Calculate the Se & Te for each region Correctly!!
    if rg == 9
        id = 1;
    elseif ( rg == 1 || rg == 2 || rg == 4 || rg == 5 )  % overlapping XY PML region
        id = 4;
    elseif ( rg == 6 || rg == 7 )                        % Y PML region
        id = 3;
    else%if ( rg == 3 || rg == 8 )                       % X PML region
        id = 2;
    end
    
    %     S = zeros(3,3);
    Te = ezz(id)*Ae * T;                     % Calculates Local Mass Matrix
    for ii = 1:3
        for jj = 1:3                         % Calculates Local Stiffness Matrix
            S(ii,jj) = (1/muyy(id) * b(ii)*b(jj) + 1/muxx(id) * c(ii)*c(jj))*Ae;
            if ( node_id(n(ii)) ~= 0 )
                if ( node_id(n(jj)) ~= 0 )
                    Aet(index(n(ii)) , index(n(jj))) = ...
                        Aet(index(n(ii)) , index(n(jj))) + S(ii,jj) - (omega^2)*Te(ii,jj);
                else
                    B(index(n(ii))) = B(index(n(ii))) - ( S(ii,jj) - (omega^2)*Te(ii,jj))*Es(n(jj) );
                end
                
            end
        end
        
    end
end
toc
X = Aet\B;                           % direct solver
for i =1:Nn
    if index(i) ~= 0                 % completing the result vector after the system solution
        Es(i) = X(index(i)) ;        % by creating a total vector including the boundaries
    end
end
E = Es+Ei ;                          % total Electric Field = Scattered and Incident
%%  FEM plots 
figure
subplot(2,2,1)                 % plot the real part of the total field
pdeplot(p/lambda,e,t,'XYdata',real(E),'Contour','off')
axis equal; axis tight;
colormap(jet);
hcb = colorbar;
hcb.Label.String = 'V/m';
xlim([x_min x_max]/lambda), ylim([y_min y_max]/lambda)  % XY-limits
xlabel('X-axis (\lambda)'), ylabel('Y-axis (\lambda)')  % XY-label
title('Real\{E\}')

subplot(2,2,2)                % plot the imaginary part of the total field
pdeplot(p/lambda,e,t,'XYdata',imag(E),'Contour','off')
axis equal; axis tight;
colormap(jet);
hcb = colorbar;
hcb.Label.String = 'V/m';
xlim([x_min x_max]/lambda), ylim([y_min y_max]/lambda)  % XY-limits
xlabel('X-axis (\lambda)'), ylabel('Y-axis (\lambda)')  % XY-label
title('Imag\{E\}')

subplot(2,2,3:4)              % plot the magnitude of the total field
pdeplot(p/lambda,e,t,'XYdata',abs(E),'Contour','off')
axis equal; axis tight;
colormap(jet);
hcb = colorbar;
hcb.Label.String = 'V/m';
xlim([x_min x_max]/lambda), ylim([y_min y_max]/lambda)  % XY-limits
xlabel('X-axis (\lambda)'), ylabel('Y-axis (\lambda)')  % XY-label
title('|E|')
hold off
sg = sgtitle({'Total Electric Field   E = E_i + E_s'});
sg.FontName    = 'Times';
%% Analytical Result
points     = 500;                                 % space discretization (points < 1000 !!!!)  <---------------------
xx         = linspace(x_min,x_max,points);        % x-coordinate
yy         = linspace(y_min,y_max,points);        % y-coordinate
tt         = 20;                                  % number of terms to be calculated in the analytical result (20 are enough) <---------------------
epsilon_nu = [1 2*ones(1,tt)];                    % coefficient needed for the calculation of the Analytical expansion
Es_Theor   = 0*zeros(length(xx),length(yy));      % scattered field

for nu = 0:tt
    temp = nan*zeros(length(xx),length(yy));                 % accoutns for nu-th term of the expansion 
    a_nu = besselj(nu,k0*Cyl_Rad)/besselh(nu,2,k0*Cyl_Rad);  % weight coefficient for every term in the expansion
    
    for ii = 1:length(yy)
        for jj = 1:length(xx)
            if sqrt(xx(jj)^2 + yy(ii)^2) >= Cyl_Rad    % The expansion has the following form
%                 rho = sqrt(xx(ii)^2 + yy(jj)^2);     % ε_n * α_n * H_n^(2)(k0*\rho) * cos(n*φ)
%                 phi = atan2(yy(ii),xx(jj));
                temp(ii,jj) =  epsilon_nu(nu+1)*((-1i)^nu)*...
                    besselh(nu,2,k0*sqrt(xx(ii)^2 + yy(jj)^2))*cos(nu*atan2(yy(ii),xx(jj)));
            end
        end
    end
    Es_Theor = Es_Theor + a_nu*temp;                   % final result. conttains all the terms 
end
Es_Theor = -Eo*Es_Theor;                               % putting the last things 
%% Result Verification 
figure
subplot(1,2,1)
pdeplot(p/lambda,e,t,'XYdata',abs(Es),'Contour','off')
axis equal; axis tight;
colormap(jet);                             %      <---------------------
hcb = colorbar;
hcb.Label.String = 'V/m';
caxis([0 inf])
xlim([x_min x_max]/lambda), ylim([y_min y_max]/lambda)  % XY-limits
xlabel('X-axis (\lambda)'), ylabel('Y-axis (\lambda)')  % XY-label
title('Numerical'),hold off

subplot(1,2,2)
surf(xx/lambda,yy/lambda,abs(Es_Theor),'EdgeColor','none');
view(2),grid off,axis equal; axis tight;                                               
colormap(jet)                               %     <---------------------
hcb = colorbar;
hcb.Label.String = 'V/m';
caxis([0 inf])                              % set the limits of the colorbar 
xlim([x_min x_max]/lambda), ylim([y_min y_max]/lambda)  % XY-limits
xlabel('X-axis (\lambda)'), ylabel('Y-axis (\lambda)')  % XY-label
title('Analytical') ,hold off
sg = sgtitle({'Result Verification  \{E_s\}'});
sg.FontName = 'Times';



