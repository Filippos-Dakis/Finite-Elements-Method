function [E2,H2,Ex,Ey,Hx,Hy] = EH_Abs_CrSection(Fx,Fy,mode,fc,freq,e0er,mu0)
%% Filippos Tzimkas-Dakis  MSc. Student, UoC Physics Dept. September 2021
%
% This function calculates the absolute value of Electric (E) and Magnetic (H) 
% field components. The calculation depends on the whether the mode is TM or TE.            
% [Fx,Fy] = gradient(Field_z_component);        | You can call the function with only these three
% mode    = 'TE' or 'TM'                        | inputs and get a "general" result
% fc      = cutoff frequency in Hz
% freq    = the frequency in which you want to see the mode profile 
% e0er    = dielectric constant in the waveguide
% mu0     = vacumm permeability

switch nargin
    case 3
        fc   = 1/(2*pi);
        freq = 1;
        e0er = 1;
        mu0  = 1;
        beta = 1;
        kc   = 1;
        w    = 1;
    case 7
        k    = 2*pi*sqrt(mu0*e0er)*freq;
        kc   = 2*pi*sqrt(mu0e0*er)*fc;
        beta = sqrt(k^2 - kc^2);
        w    = 2*pi*freq;
end

if strcmp(mode,'TE')
   % Ez = 0 
    Ex = -( 1i*w*mu0/kc^2 ) * Fy;
    Ey =  ( 1i*w*mu0/kc^2 ) * Fx;
    E2 = sqrt( abs(Ex).^2 + abs(Ey).^2 );             % calculates the field amplitude
    E2 = E2/max(max(E2));                             % normalizes amplitude with maximum value
    
    % Hz ~= 0 
    Hx =  -( 1i*beta/kc^2 ) * Fx;
    Hy =  -( 1i*beta/kc^2 ) * Fy;
    H2 = sqrt( abs(Hx).^2 + abs(Hy).^2 );             % calculates the field amplitude
    H2 = H2/max(max(H2));                             % normalizes amplitude with maximum value
elseif strcmp(mode,'TM')
    % Ez ~= 0 
    Ex =  -( 1i*beta/kc^2 ) * Fx;
    Ey =  -( 1i*beta/kc^2 ) * Fy;
    E2 = sqrt( abs(Ex).^2 + abs(Ey).^2 );             % calculates the field amplitude
    E2 = E2/max(max(E2));                             % normalizes amplitude with maximum value
    
    % Hz ~= 0 
    Hx =  ( 1i*w*e0er/kc^2 ) * Fy;
    Hy = -( 1i*w*e0er/kc^2 ) * Fx;
    H2 = sqrt( abs(Hx).^2 + abs(Hy).^2 );             % calculates the field amplitude
    H2 = H2/max(max(H2));                             % normalizes amplitude with maximum value   
else
    error('~~ ERROR:  mode input must be "TE" or "TM" !!!')
end


end

