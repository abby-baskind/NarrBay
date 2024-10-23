% o2satv2a.m                                       by:  Edward T Peltzer, MBARI
%                                                  revised:  2013 Sep 23.
%
% CALCULATE OXYGEN CONCENTRATION AT SATURATION
%
% Source:  Garcia & Gordon (1992).  Oxygen solubility in seawater:
%          Better fitting equations.  L&O 37: 1307-1312.
%
% Input:       S = Salinity (pss-78)
%              T = Temp (deg C)
%
% Output:      Oxygen saturation at one atmosphere (mL/L).
%
%                        O2 = o2satv2a(S,T).

function [O2] = o2satv2a(S,T)


% DEFINE CONSTANTS, ETC FOR SATURATION CALCULATION

%  The constants used are for units of mL O2/L.

  A0 =  2.00856;
  A1 =  3.22400;
  A2 =  3.99063;
  A3 =  4.80299;
  A4 =  9.78188e-01;
  A5 =  1.71069;
  
  B0 = -6.24097e-03;
  B1 = -6.93498e-03;
  B2 = -6.90358e-03;
  B3 = -4.29155e-03;
  
  C0 = -3.11680e-07;
  
%   Calculate Ts from T (deg C)  

  Ts = log((298.15 - T) ./ (273.15 + T));
  
%   Calculate O2 saturation in mL O2/L.

  A = ((((A5 .* Ts + A4) .* Ts + A3) .* Ts + A2) .* Ts + A1) .* Ts + A0;
  
  B = ((B3 .* Ts + B2) .* Ts + B1) .* Ts + B0;
  
  O2 = exp(A + S.*(B + S.*C0));
