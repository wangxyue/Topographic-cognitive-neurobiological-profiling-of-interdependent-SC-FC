function x = ACEcode(ACE)
%
% Maps active (non-zero) ACE parameters to a integer code
%
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/34ea883
%          2018-09-21 09:59:07 +0100

Code = [1 2 3];
x    = sum( Code(abs(ACE)>0) );

return

