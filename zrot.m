%	function [M] = zrot(angle)
%
%	Function returns the rotation matrix M such that
%	y = Mx rotates a 1x3 cartesian vector about the z axis
%	by angle degrees.
%

% ======================== CVS Log Messages ========================
% $Log: zrot.m,v $
% Revision 1.2  2002/03/28 00:50:12  bah
% Added log to source file
%
%
%
% ================================================================== 


function [M] = zrot(angle)

c = reshape(cos(pi*angle/180), 1, 1, length(angle));
s = reshape(sin(pi*angle/180), 1, 1, length(angle));

M = [c s zeros(1, 1, length(angle));
     -s c zeros(1, 1, length(angle));
     zeros(1, 1, length(angle)) zeros(1, 1, length(angle)) ones(1, 1, length(angle))];

% by Brian Hargreaves
