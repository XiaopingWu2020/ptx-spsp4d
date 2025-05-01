function tesla = hz2tesla (hz, gamma)
% HZ2TESLA Convert hz to tesla.
%
% Usage: tesla = hz2tesla (hz, gamma)
%
% Returns
% -------
% tesla: 
%
% Expects
% -------
% hz: 
% gamma: defaults to 2.675e8 rad/sec/tesla for proton.
%
%
% See also: tesla2hz
%
%
% Copyright (C) 2008 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Mon Aug 18 15:56:50 2008
%

if nargin < 2
  gamma = 2.675e8;
end

tesla = 2.*pi.* hz ./ gamma;
