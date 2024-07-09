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

function [C,h] = polarcont(theta,r,z,N,s)
[a,b] = size(z);
if a ~= length(r)
    error('r is not the same length as the first dimension of z')
end
if b ~= length(theta)
    error('theta is not the same length as the second dimension of z')
end
x = zeros(a,b);
y = zeros(a,b);
theta=90-theta;
for j = 1:a
    for k = 1:b
        x(j,k) = r(j)*cosd(theta(k));
        y(j,k) = r(j)*sind(theta(k));
    end
end
if nargin == 3
    [C,h] = contourf(x,y,z,'LineColor','none');
elseif nargin == 4
    [C,h] = contourf(x,y,z,N,'LineColor','none');
elseif nargin == 5
    [C,h] = contourf(x,y,z,N,s,'LineColor','none');
else
    error('Incorrect number of inputs')
end