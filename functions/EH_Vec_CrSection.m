function [Ex,Ey,Hx,Hy] = EH_Vec_CrSection(p,e,t,dl,units,X,mode,fc,freq,e0er,mu0)
%% Filippos Tzimkas-Dakis  MSc. Student, UoC Physics Dept. September 2021
%
% This fuction servers as a tool to calculate and plot the dynamic lines of
% the Electric and Magnetic field inside a metallic waveguide. This function
% calculates all the E/H components on a cross-section and then plots the
% results. One plot for [Ex,Ey] and one for [Hx,Hy].
% Attention!! You cannot use these components to calculate the magnitude of
% the fields!! You will get a wrong answer !!
%
% p,e,t,dl have to do with the geometry of the waveguide
% units = 1 --> meters,  0.01 --> cm ,  0.001 --> mm etc
% X     = the Ez(TM-mode) or Hz(TE-mode) components calculated through FEM
% mode  = either 'TE' or 'TM'
% fc    = the cutoff frequency calculated from the eigenvalue problem
% freq  = the frequecy we are interested in
% e0er  = the dielectric constant of the medium that fills the waveguide
% mu0   = the vacumm permeability
%
% For more detail read  Chapter 3 from "Microwave Engineering" (David Pozar)
%
% Do not hesitate to conatct me at    
%           dakisfilippos@gmail.com    or   f.tzimkas@physics.uoc.gr
%
switch nargin
    case 7
        fc   = 1/(2*pi);
        freq = 1;
        e0er = 1;
        mu0  = 1;
        beta = 1;
        kc   = 1;
        w    = 1;
    case 11
        k    = 2*pi*sqrt(mu0*e0er)*freq;
        kc   = 2*pi*sqrt(mu0e0*er)*fc;
        beta = sqrt(k^2 - kc^2);
        w    = 2*pi*freq;
end
[Fx,Fy] = pdegrad(p,t,X);
p     = p/units;
x_max = max(p(1,:));                 % max x-position
x_min = min(p(1,:));                 % min x-position
y_max = max(p(2,:));                 % max y-position
y_min = min(p(2,:));                 % min y-position
units_c = [1, 1e-1 1e-2 1e-3 1e-6 1e-9];       % setting the graph units 
units_C = {'m','dm','cm','mm','um','nm'};
index    = find(units_c == units);
un     = units_C{index};
x_axis = ['X-axis (',un,')'];
y_axis = ['Y-axis (',un,')'];
figure
if strcmp(mode,'TE')
    % Ez = 0
    Ex = -( 1*w*mu0/kc^2 ) * Fy;
    Ey =  ( 1*w*mu0/kc^2 ) * Fx;
    subplot(1,2,1)
    pltE = pdeplot(p,e,t,'FlowData',[Ex;Ey]);          % plots the Electric Field
    pltE.Color =  [1 1 0];
    hold on
    set(gca,'Color','k')                               % sets the background color black
    hold on
    h = pdegplot(dl);                                  % plots the resulted geometry
    h.LineWidth = 1.5;
    h.Color = [1 1 1];                                 % with white edges
    box on
    axis equal
    title('Electric Field')
    xlabel(x_axis)
    ylabel(y_axis)
    xlim([x_min x_max]*1.1)
    ylim([y_min y_max]*1.1)
    hold off
    
    subplot(1,2,2)
    % Hz ~= 0
    Hx =  ( 1*beta/kc^2 ) * Fx;
    Hy =  ( 1*beta/kc^2 ) * Fy;
    pltE = pdeplot(p,e,t,'FlowData',[Hx;Hy]);          % plots the Electric Field
    pltE.Color = [1 1 0];
    hold on
    set(gca,'Color','k')                               % sets the background color black
    hold on
    h = pdegplot(dl);                                  % plots the resulted geometry
    h.LineWidth = 1.5;
    h.Color = [1 1 1];                                 % with white edges
    box on
    axis equal
    title('Magnetic Field')
    xlabel(x_axis)
    ylabel(y_axis)
    xlim([x_min x_max]*1.1)
    ylim([y_min y_max]*1.1)
    hold off
    
elseif strcmp(mode,'TM')
    % Ez ~= 0
    Ex =  -( 1*beta/kc^2 ) * Fx;
    Ey =  -( 1*beta/kc^2 ) * Fy;
    subplot(1,2,1)
    pltE = pdeplot(p,e,t,'FlowData',[Ex;Ey]);          % plots the Electric Field
    pltE.Color =  [1 1 0];
    hold on
    set(gca,'Color','k')                               % sets the background color black
    hold on
    h = pdegplot(dl);                                  % plots the resulted geometry
    h.LineWidth = 1.5;
    h.Color = [1 1 1];                                 % with white edges
    box on
    axis equal
    title('Electric Field')
    xlabel(x_axis)
    ylabel(y_axis)
    xlim([x_min x_max]*1.1)
    ylim([y_min y_max]*1.1)
    hold off
    
    subplot 122
    % Hz = 0
    Hx = -( 1*w*e0er/kc^2 ) * Fy;
    Hy =  ( 1*w*e0er/kc^2 ) * Fx;
    subplot(1,2,2)
    pltE = pdeplot(p,e,t,'FlowData',[Hx;Hy]);          % plots the Electric Field
    pltE.Color = [1 1 0];
    hold on
    set(gca,'Color','k')                               % sets the background color black
    hold on
    h = pdegplot(dl);                                  % plots the resulted geometry
    h.LineWidth = 1.5;
    h.Color = [1 1 1];                                 % with white edges
    box on
    axis equal
    title('Magnetic Field')
    xlabel(x_axis)
    ylabel(y_axis)
    xlim([x_min x_max]*1.1)
    ylim([y_min y_max]*1.1)
    hold off
else
    error('~~ ERROR:  mode input must be "TE" or "TM" !!!')
end
sg = sgtitle({[mode,' mode'],'\ $\vec{\rm E}_t$ \& $\vec{\rm H}_t$ profiles'},'Interpreter','latex');
sg.FontName    = 'Times';
end

