function Fbub = bubbleinjection(wind,T,DO,DOsat)

% calculate bubble injection based on Bushinsky and Emerson, 2018
% Biological and physical controls on the oxygen cycle in the 
% Kuroshio Extension from an array of profiling floats
% 
% F_tot = Fs + Fc + Fp
% 
% Fs: flux from air-sea interface exchange
% Fc: flux from small collapsing bubbles
% Fp: flux from large, partially collapsing bubbles
% 
% DO and DOsat unit: mg/L = g/m3
% unit transfer to mol/m3 in calculation
% 
% See Section 4.1.1 in Bushinsky & Emerson (2018)
% REFERENCE: Bushinsky, S. M., & Emerson, S. R. (2018). 
% % Biological and physical controls on the oxygen cycle in the Kuroshio 
% % Extension from an array of profiling floats. 
% % Deep Sea Research Part I: Oceanographic Research Papers, 141, 51-70.
% 
% Written by Shuai Gu, February 2024.
% Edited by Abby Baskind, March 2024.


% DRAG COEFFICIENT *******************************************************
% wind speed [m/s]    ||    Cd [unitless]
% 0-11                ||    0.0012
% 11-20               ||    (0.49 + 0.065 * U10) * 1e-3
% 18+                 ||    0.0018
Cd=[];

for i=1:length(wind)
   
    if wind(i)<11;
        Cd(i,1) = 0.0012;
    elseif wind(i)>=11 && wind(i)<=20;
        Cd(i,1) = (0.49+0.065.*wind(i)).*10^-3;
    elseif wind(i)>20;
        Cd(i,1) = 0.0018;
    elseif isnan(wind(i));
        Cd(i,1)=nan;
    end
    
end

% SCHMIDT NUMBER *********************************************************
% Schmidt Number (Wanninkhof, 1992)
Sc =(1953.4 - 128 .* T + 3.9918 .* T.^2 - 0.050091 .*T.^3);

% AIR SIDE FRICTION VELOCITY *********************************************
% Ua_star = Cd^(0.5) * U10 [m/s]
u_star_a = wind .* sqrt(Cd);                                % unit: m/s

% Water SIDE FRICTION VELOCITY *******************************************
% Uw_star = (rho_a/rho_w)^(0.5) * Ua_star [m/s]
u_star_w = 0.034 .* u_star_a;                               % unit: m/s

% FS *********************************************************************
% Fs = ks(O2sat - O2)
% ks: air sea interface mass transfer coeff
% % ks = 1.3e-4 * Ua_star * (Sc/660)^(-0.5)
ks = 0.00013 .* u_star_a .* ((Sc/660).^(-0.5));             % unit: m/s
ks_factor = 1;
Fs = ks .* ks_factor .* (DOsat - DO) .* 86400;              % unit gO2/m2d 

% FC: SMALL BUBBLE FLUX **************************************************
% "Small bubbles that completely collapse add gas to the water 
% as a function of wind speed but are independent of saturation state."
% Bushinsky & Emerson
% Fc = 5.56 * B * Uw_star^(3.86) * X_O2
% % X_O2 = 0.20946: atmosphere mole fraction of oxygen
kc = 5.56 .* (u_star_w.^3.86); 
kc_factor = 1;
Fc = kc_factor .* kc .* 0.20946 * 86400 * 32;               % unit: gO2/m2d

% FP: LARGE BUBBLE FLUX **************************************************
% Fp = 5.5B * Uw_star^(2.76) * (Sc/660)^(-2/3) * "stuff"
% "stuff" = [1 + ∆p] * O2sat - O2
% ∆p = 1.5244 * Uw_star^(1.06): fractional increase in pressure
% experienced by bubbles driven below the surface due to breaking 
% waves at high wind speeds

d_P = 1.5244 .* (u_star_w.^1.06);
kp_factor = 1;
kp = 5.5 .* (u_star_w.^2.76) .* ((Sc/660).^(-2/3));
Fp = kp_factor .* kp .* ((1+d_P) .* DOsat - DO) .* 86400;   % g m-2 day-1

% TOTAL BUBBLE INJECTION *************************************************
% Fc + Fp = little bubbles + big bubbles
% unit: g m-2 d-1
Fb = Fc + Fp;

Fbub = Fb;
