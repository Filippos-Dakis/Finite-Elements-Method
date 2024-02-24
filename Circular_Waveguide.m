%% Filippos Tzimkas-Dakis  MSc. Student, UoC Physics Dept. September 2021
%
% I wrote this short because I wanted to learn some things on the FEM
% computational method. There is a Circular waveguide (infinite length),
% with radius R. I solve the eigenvalue problem to obtain the cutoff
% frequencies and also modes profile, both for TE and TM modes.
% The user can chnage several things, such as the number of refinements,
% the number of modes to be calculated, the mode he/she wants to focus on
% (calculate all the E/M components ) and many more.
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
close all,clear all, clc
%                                     choose the mode you are interested in calcualting
MODE = 'TM';                        % 'TE' or 'TM'                                     <---------------------
N_modes = 4;                        % choose the number of modes you want to calculate <---------------------
cm    = 1e-2;                       % set the length units of the radii you will use   
mm    = 1e-3;                       % m->1, cm->1e-2, mm-> 1e-3 etc
units = cm;                         % change the units of the graphs !!                <---------------------

R_out = 1;                          % (m_units) outer radius                           <---------------------
R_out = R_out * cm;                 % setting the dimensions correctly                 <---------------------
e0    = 8.8541878128*1e-12;         % (F/m) Vacuum permitivity
e_r   = 1;                          % (clean) relative permitivity                     <---------------------
mu0   = 4*pi*1e-7;                  % (H/m) vacuum permiability

% creating the geometry. search for 'decsg' in MATLAB documentation
C1 = [1 0 0 R_out]';                % outter circle; column vector
gd = [C1];                          % needed for matlab
ns = char('C1');                    % names of the boundaries
ns = ns';                           % converting it to column vector
sf = 'C1';                          % setting the formula for the geometry

[dl,bt] = decsg(gd,sf,ns);          % creates the geometry

% pdegplot(dl,'EdgeLabels','on','FaceLabels','on');   % plots the resulted geometry
% xlim([-R_out,R_out])
% axis equal
%% plot parameters for later on
units_c = [1, 1e-1 1e-2 1e-3 1e-6 1e-9];       % setting the graph units
units_C = {'m','dm','cm','mm','um','nm'};
index   = find(units_c == units);
un      = units_C{index};
x_axis  = ['X-axis (',un,')'];                 % labels for plots later on
y_axis  = ['Y-axis (',un,')'];                 % labels for plots later on
%% -------------- creating the mesh  --------------------------------------
refin   = 1;                         % choose the number of refinments   <--------------------------------------

[p,e,t] = initmesh(dl);              % creates the triangularization
for ii = 1:refin
    [p,e,t] = refinemesh(dl,p,e,t);  % denses the triangular lattice
end
pdeplot(p/units,e,t);                % plots the geometry with the triagualrization
xlabel(x_axis),ylabel(y_axis)
title(['Lattice after ',num2str(refin),' refinement(s)']), hold off
%% 2 "Re-number" the unknown nodes ----------------------------------------
Nn      = size(p,2);                 % number of nodes
Ne      = size(t,2);                 % number of elements
Nd      = size(e,2);                 % number of edges - ακμές
node_id = ones(Nn,1);                % our flag in order to check if the node
%                                      is on the edge or not
X0      = zeros(Nn,1);               % initialization of the result vector
% In this part we apply the boundary conditions.
if strcmp(MODE,'TM')                 % This If checks if we have the user
    Bound_Cond = 1;                  % has choosen TE or TM mode, if not the program stops.
    profls     = 'E_z';
    fprintf('\n');
    cprintf('*yellow',' TM Mode Operation \n')
    cprintf('*yellow',' Dirichlet Boundary Conditions \n')
elseif strcmp(MODE,'TE')             % TE mode
    Bound_Cond = 0;
    profls     = 'H_z';
    fprintf('\n');
    cprintf('*cyan',' TE Mode Operation \n')
    cprintf('*cyan',' von Neumann Boundary Conditions \n')
else                                 % stops if MODE ~= TE  || TM
    error('~~ QUACK! You have to chose MODE = "TE" or "TM" !! ~~')
end

for id = 1:Nd
    if ( e(6,id )== 0 || e(7,id) == 0 )
        if Bound_Cond
            node_id(e(1,id)) = 0;                % marks the specific nodes as KNOWN because they are on the edges
            node_id(e(2,id)) = 0;                % end we have to apply the boudnary conditions
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
%% 3. Main - Kernel Calculation of Stiffness and Mass Matrices ------------
Nf  = nnz(node_id);                 % number of uknown variables
Tff = spalloc(Nf,Nf,7*Nf);          % initialization of the Large Mass Matrix
Sff = spalloc(Nf,Nf,7*Nf);          % initialization of the large Stiff Matrix for system Solutions Sff*X = B
Te  = zeros(3,3);                   % initialization of the Local Mass Matrix
Se  = zeros(3,3);                   % initialization of the Local Stiff Matrix

Te = (1/12) * ones(3);                  % local Mass Matrix is standard !
Te = Te - (1/12)*eye(3) + (1/6)*eye(3); % Te = Ae* [1/6    1/12   1/12;
%                                                  1/12   1/6    1/12;
%                                                  1/12   1/12    1/6 ];

for ie = 1:Ne              % Scans all the Elements (triangles)
    n(1:3) = t(1:3,ie);            % n(1:3) are the nodes of every elemnet
    rg     = t(4,ie);              % rg is the region of the element (subdomain number)
    x(1:3) = p(1,n(1:3));          % stores the x-coordinate
    y(1:3) = p(2,n(1:3));          % stores the y-coordinate
    De     = det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3)]);  % determinant
    Ae     = abs(De/2);            % element area
    
    b = zeros(3,1);
    c = zeros(3,1);
    for k = 0:2
        xk = circshift(x,-k);         % Coefficients for linear
        yk = circshift(y,-k);         % triangular finite elements. For more
        b(k+1) = (yk(2)-yk(3))/De;    % information on these equations read Tsimboukis'
        c(k+1) = (xk(3)-xk(2))/De;    % complementary books "Notes on Computational Electromagnetics" pages 83-95
    end                               % and "Energy Methods for E/M Fields" pages 287-308
    
    for ii = 1:3                      % calculates the 3X3 Stiff Local Matrix
        for j = 1:3                   %
            Se(ii,j) = ( b(ii)*b(j) + c(ii)*c(j) )*Ae;  % Local Stiff Matrix
            if ( node_id(n(ii)) ~= 0 )                      % node n(i) is Unknown
                if ( node_id(n(j)) ~= 0 )                   % node n(j) is Unknown
                    Sff(node_id(n(ii)),node_id(n(j))) = Sff(node_id(n(ii)),node_id(n(j))) + Se(ii,j);
                    Tff(node_id(n(ii)),node_id(n(j))) = Tff(node_id(n(ii)),node_id(n(j))) + Te(ii,j)*Ae;
                end
            end
        end
    end
end
%% 4. Calculates the Eigenvalues and Eigenvectors for the system (S - T)x = 0
[V,D] = eigs(Sff,Tff,N_modes,'smallestabs');   % Solves the eigenvalue problem
D     = diag(D);                               % creates a vector
fc    = sqrt(D)/(2*pi*sqrt(e_r*e0*mu0));       % Cutoff frequencies
pnm   = sqrt(D)*R_out                          % prints the eigenvalues in command window   
%% Converts the geometry into "units"meter   (m,cm,mm etc)
dl_plot = dl;
dl_plot(2:end,1:4) = dl(2:end,1:4)/units;
%% 5. Plots modes profiles TE(E_t = E_x + E_y) and TM(H_t = H_x + H_y)
fig1          = figure(2);            % To make sure that the saved PDF will be 100% vector      <------ These commands solve the problem of
fig1.Renderer = 'Painters';           % because I hade some problems !!

if sqrt(N_modes) == fix(sqrt(N_modes))
    rowz = fix(sqrt(N_modes));
    colmz = fix(N_modes/rowz);
else
    rowz = fix(sqrt(N_modes)) + 1;
    colmz = fix(N_modes/rowz) + 1;
end

X1 = zeros(Nn,length(D));
for ii = 1:length(fc)
    X = V(:,ii);
    % complete the X0 matrix with the computed values of the unknown nodes
    for k = 1:Nn
        if( node_id(k) > 0 )
            X0(k) = X(node_id(k));
        end
    end
    X1(:,ii) = X0;
    subplot(rowz,colmz,ii)
    pltF = pdeplot(p/units,e,t,'XYdata',X0,'Contour','on');  % plots the Potential in colormap
    set(gca,'Color','k'),hold on
    h = pdegplot(dl_plot);                                  % plots the resulted geometry
    h.Color = [1 1 1];                                      % with white color
    title(['f_c = ', num2str(round(fc(ii)*1e-9,2)), 'GHz'])
    colormap(viridis);
    axis equal, axis tight, box on
end
sg = sgtitle({[MODE,' Modes, ',profls,'(x,y)'],['_{( x,y axes in [',un,'] )}']});
sg.FontName = 'Times';

Field_z = X1;               % store the resulting Z components (Ez for TM and Hz for TE)
fc_TM   = fc*1e-9;          % cutoff frequencies in GHz
%% Plot the profile (Magnitude and Vector) of a choosen Mode
%  Check the figure which contains all the calculated modes and chose the
%  one you are interested and plote [Ex,Ey] and [Hx,Hy]
E_H     = X1(:,2);         % choose the mode you want to see in detail  <---------------------

F       = pdeInterpolant(p,t,E_H);                    % preparation for interpolation command
XX1     = linspace(-R_out,R_out,1000);                % points along X-axis
YY1     = linspace(-R_out,R_out,1000);                % points along Y-axis
dX1     = abs(XX1(1) - XX1(2));                       % x-step
dY1     = abs(YY1(1) - YY1(2));                       % y-step
[Xx1,Yy1] = meshgrid(XX1,YY1);                        % grid generation
uOut    = evaluate(F,Xx1,Yy1);                        % interpolation command execution
Z       = reshape(uOut,size(Xx1,1),size(Yy1,1));      % from vector to array
[Fx,Fy] = gradient(Z);                                % calcultes the Electric Field
Fx      = Fx/dX1;                                     % inserting the x-step to derivative
Fy      = Fy/dY1;                                     % inserting the y-step to derivative

Z2      = sqrt( abs(Fx).^2 + abs(Fy).^2 );            % calculates the field amplitude
Z2      = Z2/max(max(Z2));

[E2,H2,Ex,Ey,Hx,Hy] = EH_Abs_CrSection(Fx,Fy,MODE);

x_max = max(p(1,:));                 % max x-position
x_min = min(p(1,:));                 % min x-position
y_max = max(p(2,:));                 % max y-position
y_min = min(p(2,:));                 % min y-position

figure
subplot(3,2,1)
surf(Xx1/units,Yy1/units,E2,'EdgeColor','none')
view(2)                                              % 2D view
colormap(viridis),colorbar
set(gca,'Color','k'),hold on                         % sets the background Black
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                   % with white color
xlim([x_min x_max]*1.1/units),ylim([y_min y_max]*1.1/units) % XY-limits
xlabel(x_axis),ylabel(y_axis)
title('Electric Field')

subplot(3,2,3)
surf(Xx1/units,Yy1/units,imag(Ex)/max(max(imag(Ex))),'EdgeColor','none')
view(2)                                              % 2D view
colormap(viridis),colorbar
set(gca,'Color','k'),hold on                         % sets the background Black
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                   % with white color
xlim([x_min x_max]*1.1/units),ylim([y_min y_max]*1.1/units) % XY-limits
xlabel(x_axis),ylabel(y_axis)
title('Im\{E_x\}')

subplot(3,2,5)
surf(Xx1/units,Yy1/units,imag(Ey)/max(max(imag(Ey))),'EdgeColor','none')
view(2)                                              % 2D view
colormap(viridis),colorbar
set(gca,'Color','k'),hold on                         % sets the background Black
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                   % with white color
xlim([x_min x_max]*1.1/units),ylim([y_min y_max]*1.1/units) % XY-limits
xlabel(x_axis),ylabel(y_axis)
title('Im\{E_y\}')

subplot(3,2,2)
surf(Xx1/units,Yy1/units,H2,'EdgeColor','none')
view(2)                                              % 2D view
colormap(viridis),colorbar
set(gca,'Color','k'),hold on                         % sets the background Black
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                   % with white color
xlim([x_min x_max]*1.1/units),ylim([y_min y_max]*1.1/units) % XY-limits
xlim([x_min x_max]*1.1/units),ylim([y_min y_max]*1.1/units) % XY-limits
xlabel(x_axis),ylabel(y_axis)
title('Magnetic Field')

subplot(3,2,4)
surf(Xx1/units,Yy1/units,imag(Hx)/max(max(imag(Hx))),'EdgeColor','none')
view(2)                                              % 2D view
colormap(viridis),colorbar
set(gca,'Color','k'),hold on                         % sets the background Black
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                   % with white color
xlim([x_min x_max]*1.1/units),ylim([y_min y_max]*1.1/units) % XY-limits
xlabel(x_axis),ylabel(y_axis)
title('Im\{H_x\}')

subplot(3,2,6)
surf(Xx1/units,Yy1/units,imag(Hy)/max(max(imag(Hy))),'EdgeColor','none')
view(2)                                              % 2D view
colormap(viridis),colorbar
set(gca,'Color','k'),hold on                         % sets the background Black
h = pdegplot(dl_plot);                               % plots the resulted geometry
h.LineWidth = 1.5;
h.Color = [1 1 1];                                   % with white color
xlim([x_min x_max]*1.1/units),ylim([y_min y_max]*1.1/units) % XY-limits
xlabel(x_axis),ylabel(y_axis)
title('Im\{H_y\}')
    
sg = sgtitle({[MODE,' mode'],' $\rm {|\vec{E}_t|}_{normalized}$ \& $\rm {|\vec{H}_t|}_{normalized}$'},'Interpreter','latex');
sg.FontName    = 'Times';

% Vector plots at the cross-section
% read the documentation of EH_Vec_CrSection(p,e,t,dl,X,mode,fc,freq,e0er,mu0)
[Ex,Ey,Hx,Hy] = EH_Vec_CrSection(p,e,t,dl_plot,units,E_H,MODE);


