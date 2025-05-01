function sysmat = construct_sysmat_spspkT (kp, b1maps, mask, fox, b0map, ...
                                           phasetrack, freqs, wts, dt, poffset)

% CONSTRUCT_SYSMAT_SPSPKT Construct system matrix for kT point based spsp pulse design.
%
% Usage: sysmat = construct_sysmat_spspkT (kp, b1maps, mask, fox, b0map,
% phasetrack, freqs, wts, dt, poffset)
%
% Returns
% -------
% sysmat: output system matrix.
%
% Expects
% -------
% kp: 3-by-nktpts matrix for kT points in kspace, rad/m
% b1maps: b1 maps
% mask: spatial mask
% fox: field of excitation in m.
% b0map: b0 map, will be zero when its empty or not specified.
% 
% phasetrack: a vector of time points indicating when to trace the phase due to
% off resonance (i.e., b0 and/or chemical shifts). The number of its elements
% should equal the number of the kT points or of the subpulses. defaults to
% 1:nktpts.
% 
% freqs: a vector of chemical shift related freq offsets in hz that are of
% interest.
% 
% wts: a vector of weights used to control the relative importance of freqs.
% defaults to a unity vector, i.e., ones(size(freqs)).
% 
% dt: dwell in s. defaults to 10e-6.
% 
% poffset: [offsetx,offsety,offsetz] in mm specifying the offset of FOV with
% respect to grad isocenter. defaults to [0 0 0].
%
%
% See also: construct_sysmat_kTpoint construct_sysmat_spsp3d
%
%
% Copyright (C) 2011 CMRR at UMN
% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
% Created: Mon Nov 21 10:40:03 2011
%

if nargin< 10
  poffset = [0 0 0];
end
if nargin< 9
  dt = 10e-6;
end
if nargin< 8|| isempty(wts)
  wts = ones(size(freqs));
end
if nargin < 5|| isempty(b0map)
  b0map = zeros(size(mask));
end
if nargin < 6|| isempty(phasetrack)
  phasetrack = 1:size(kp,2);
end

nfreqs = length(freqs);
nchs = size(b1maps,4);
nt = size(kp,2);
nspts = length(mask(mask));

disp('-> Constructing system matrix...')
sysmat = complex(zeros(nfreqs*nspts,nchs*nt));
for ifreq= 1:nfreqs,
  iIdx0 = (ifreq-1)*nspts + 1;
  sysmat(iIdx0:iIdx0+nspts-1,:) = wts(ifreq)*constructBaseSysmat(kp, b1maps, ...
                                                    mask, fox, b0map+ ...
                                                    hz2tesla(freqs(ifreq)), ...
                                                    phasetrack,dt,1e-3* ...
                                                    poffset);
end

disp('-> System matrix created...')

%% ====================
%%
%%  Local functions
%%
%% ====================

function sysmat = constructBaseSysmat(kp,b1maps,mask,fox,offres,phasetrack, ...
                                      dt,poffset)
% 

gamma = 2.675e8;

[b1arr,posarr] = create_array(b1maps,mask,fox,poffset);

gObj = gradPulse(ones(1,phasetrack(end)),dt,'unitgrad','Tx');
kb0 = gObj.KspaceTrajectory;

kernalmat= 1i.*gamma.*dt.*exp(1i*(posarr*kp + offres(mask)*kb0(phasetrack) ) );

nchs = size(b1arr,2);
[nspa,nspo] = size(kernalmat);
sysmat = complex(zeros(nspa,nchs*nspo)); % this prealloc results in much faster construct
for idx = 1:nchs,
  iidx0 = (idx-1)*nspo + 1;
  sysmat(:,iidx0:(iidx0+nspo-1)) = diag(b1arr(:,idx)) * kernalmat;
end
