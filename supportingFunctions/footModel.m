function z = footModel( phi, legs, order )
%FOOTMODEL generate Blaberus discoidalis foor trajectories based on Fourier
% series model fitted to experimental data. 
%  z = footModel( phi, legs, order )
%    phi - phase in cycle (radians). Size Nx1 or NxM if only one leg used
%    legs - which legs to generate data for (default 1:6)
%        legs are FMH Right, HMF Left as 1-6
%    order - model order (1..4), negative values give velocities (default 3)
%    z - positions (velocities) in units of bodylength (per radian)
%
% $Revision: 1.1 $
% PolyPEDAL Lab, Berkeley 2006 by Shai Revzen (funded by NSF FIBR)
coef = [ ...
   0.4560 - 0.2480i   0.1350 - 0.3870i  -0.3160 - 0.3130i  -0.3160 + 0.3130i   0.1350 + 0.3870i   0.4560 + 0.2480i
   0.0620 + 0.0180i  -0.0720 + 0.0080i   0.0540 - 0.0340i  -0.0670 + 0.0030i   0.0670 + 0.0010i  -0.0580 - 0.0120i
   0.0600 - 0.0130i  -0.0670             0.0660 + 0.0020i  -0.0550 - 0.0330i   0.0730 + 0.0090i  -0.0600 + 0.0190i
  -0.0080 + 0.0030i   0.0030 + 0.0190i   0.0090 - 0.0010i  -0.0080 + 0.0090i  -0.0050 + 0.0120i  -0.0050 + 0.0060i
  -0.0040 - 0.0060i  -0.0060 - 0.0120i  -0.0080 - 0.0100i   0.0090 + 0.0010i   0.0040 - 0.0190i  -0.0090 - 0.0030i
   0.0070 - 0.0030i   0.0020             0.0020 + 0.0030i        0 + 0.0070i   0.0030 - 0.0020i  -0.0030          
   0.0040 - 0.0010i  -0.0030 - 0.0040i  -0.0020 + 0.0050i  -0.0030 + 0.0040i  -0.0020 + 0.0020i  -0.0050 - 0.0020i
   0.0020 + 0.0020i   0.0020            -0.0020 - 0.0010i   0.0010             0.0010 + 0.0030i        0 + 0.0010i
        0 - 0.0010i        0 - 0.0030i   0.0010            -0.0020 + 0.0010i   0.0020             0.0020 - 0.0020i
 ]';
% Coefficients are slightly different when I run it here. Ask Shai about
% the correspondence of the phi variable to the stance and swing phases.
 if nargin<3
     order = 3;
 elseif abs(order)>floor(size(coef,2)/2)
     error 'Maximal model order exceeded';
 end
 if nargin<2
     legs=[1:6];
 end
 legs = reshape(legs,1,numel(legs)); m=length(legs);
 sz = size(phi);
 if length(sz)>2
     error 'Phase values must be a vector, a matrix or a scalar'
 end
 phi = reshape(phi,1,numel(phi)); % n=length(phi);
 
 om = ([-i;i]*[0:abs(order)]);
 if order>0 
     z = coef(legs,1:2*order+1) * exp(om(2:end).'*phi);
 else
     z = coef(legs,1:-2*order+1)*diag(om(2:end))*exp(om(2:end).'*phi)./(2*pi);
 end
 if min(sz)>1
     z = reshape( z, [sz m] );
 end

% Previous model coefficients:
%    0.4780 - 0.2540i   0.1490 - 0.3950i  -0.3110 - 0.3200i  -0.3110 + 0.3200i   0.1490 + 0.3950i   0.4780 + 0.2540i
%    0.0620 + 0.0200i  -0.0740 + 0.0060i   0.0550 - 0.0340i  -0.0690 + 0.0020i   0.0690 + 0.0020i  -0.0590 - 0.0140i
%    0.0610 - 0.0150i  -0.0680 + 0.0010i   0.0680 + 0.0010i  -0.0570 - 0.0330i   0.0750 + 0.0080i  -0.0600 + 0.0210i
%   -0.0090 + 0.0030i   0.0020 + 0.0200i   0.0090 - 0.0010i  -0.0090 + 0.0090i  -0.0060 + 0.0120i  -0.0050 + 0.0060i
%   -0.0050 - 0.0060i  -0.0060 - 0.0120i  -0.0090 - 0.0090i   0.0080 + 0.0010i   0.0030 - 0.0200i  -0.0090 - 0.0030i
%    0.0070 - 0.0020i   0.0020             0.0010 + 0.0030i        0 + 0.0070i   0.0030 - 0.0020i  -0.0020 - 0.0010i
%    0.0040 - 0.0020i  -0.0030 - 0.0040i  -0.0010 + 0.0060i  -0.0030 + 0.0040i  -0.0020 + 0.0020i  -0.0050 - 0.0010i
%    0.0020 + 0.0020i   0.0020 + 0.0010i  -0.0020 - 0.0010i   0.0010 + 0.0010i   0.0010 + 0.0040i        0 + 0.0010i
%         0 - 0.0010i        0 - 0.0040i   0.0010 - 0.0010i  -0.0020 + 0.0010i   0.0020             0.0010 - 0.0020i
 
 