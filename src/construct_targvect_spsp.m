function targvect = construct_targvect_spsp (targetPat, mask, spect)

% CONSTRUCT_TARGVECT_SPSP Construct the target vector for spsp pulse design.
%
% Usage: targvect = construct_targvect_spsp (targetPat, mask, spect)
%
% Returns
% -------
% targvect: target vector assuming 1 deg flip angles.
%
% Expects
% -------
% targetPat: base target pattern
% mask: spatial mask
% 
% spect: spectrum for the bandwidth of interest.
%
%
% See also: construct_target_vector
%
%
% Copyright (C) 2011 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Mon Nov 21 11:59:05 2011
%

nfreqs = length(spect);

m = construct_target_vector(targetPat,mask,1,'s');
m = m(:);
npts = length(m);
targvect = zeros(nfreqs* npts,1);
for ifreq=1:nfreqs,
  iIdx0 = (ifreq-1)*npts + 1;
  targvect(iIdx0:iIdx0+npts-1) = spect(ifreq)* m;
end

disp('-> done...') 
