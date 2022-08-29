function [deltay, ddeltayddeltam] = UTS_CRIT(pars, Rebar, dRebar_ddeltam)
% Function for calculating the critical interference that causes the
% maximum stress to exceed the UTS strength of the material. 
%
% Inputs:
%   pars.UTS - ultimate tensile strength of material (Pa)
%   Rbar - radius after flattening of the asperity
%   dRbarddeltam - derivative of the radius with respect to maximum
%                   displacement. 

    C = 1.295*exp(0.736*pars.nu);

    deltay = (pi * C * pars.UTS/2/pars.Estar)^2*Rebar;
    
    ddeltayddeltam = (pi * C * pars.UTS/2/pars.Estar)^2*dRebar_ddeltam;

end