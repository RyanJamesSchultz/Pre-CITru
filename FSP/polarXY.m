% POLARCONT Polar contour plot
%
% Richard Rieber
% rrieber@gmail.com
% April 4, 2007
% Updated June 15, 2007
% 
% function [C,h] = polarcont(r,theta,z,N,s)
%
% Purpose: This function creates polar contour plots on the current active
%          figure
% 
% Inputs:   r     - Radius vector of length m
%           theta - Angle vector in degrees of length n
%           z     - Magnitude at the points specified in r and theta of
%                    size m x n
%           N     - The number of contours to plot [OPTIONAL]
%           s     - Linespec as described in PLOT [OPTIONAL]
%
% Outputs:  C     - returns contour matrix C as described in CONTOURC
%           h     - Column vector H of handles to LINE or PATCH objects,
%                    one handle per line.  
%
% OTHER NOTES:
% - Both C and h can be used as inputs to CLABEL
% - Colors are defined in colormap
% - Treat this function as a standard contour plot

function [h] = polarXY(theta,r,s,N,V)
theta=90-theta;
x = r.*cosd(theta);
y = r.*sind(theta);

if nargin == 2
    h = plot(x,y);
elseif nargin == 3
    h = plot(x,y,s);
elseif nargin == 5
    h = plot(x,y,s,N,V);
else
    error('Incorrect number of inputs')
end