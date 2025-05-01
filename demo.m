%%% this is a script demonstrating how to design a 4D (1D spectral and 3D spatial) 
%%% kT point pTx pulse
%%% that can be used to create a uniform water excitation (w/o touching fat spins) 
%%% across a region of interest (e.g., the entire brain). 
%%% The pTx pulse is designed in the small tip angle regime by solving a
%%% magnitude least squares minimization problem. 
%%% For a prescribed number of groups of kT points to be used, 
%%% the placement of the kT points in the excitation kspace is optimized to some extent 
%%% using a heuristic method. 
%%% 
%%% An abstract demonstrating the utility of such ptx pulses for uniform water excitation 
%%% across the brain at 7T was presented at 2012 ISMRM workshop on fat
%%% water imaging. 

clearvars; 
close all

addpath('./src/');
addpath('./regtools/')

%
load pTxCalibData.mat;%fmapMS % load 3D B1+ and B0 mapping and brain masking
mask= pTxCalibData.mask;
b1maps= 1e-6* pTxCalibData.b1maps_uTperV;
b0map= 1e-6* pTxCalibData.b0map_uT;
sliceLocations= pTxCalibData.sliceLocations_mm;
fov= pTxCalibData.FOV_mm;
res= pTxCalibData.sliceThickness_mm;
x0= pTxCalibData.phases4flipAngleMapping;

%bDrawROI= false; % true for manually drawing brain masks in slices of interest. 
%%

B0= 10.5; % in Tesla
gamma= 42.58; % Mhz/Tesla;
csfat= 3.5; % fat chemical shift in ppm. 
freqFat= round(gamma* B0* csfat);

fa= 10; % nominal flip angle in deg. 

df= 100; % frequency steps in Hz used to prescribe the frequency response.
nfp= 3;%5; % number of frequency points in water and fat bands. 
% design spectral select
waterbnd = -0.5*(df*(nfp-1)):df: 0.5*(df*(nfp-1)); %-250:125:250;
fatbnd = waterbnd- freqFat;
freqs = [fatbnd waterbnd];
wsb = 1;
wpb = 1;

%
[~,~,kdx]= ind2sub(size(mask), find(mask)); % used to determine the range of mask in the slice direction

dt = 10e-6; % s
ns = 7;%9; %size(maskMS,3); % number of slices for pulse design. 
nslicesPerSlab= (max(kdx)- min(kdx)+ 1)./ ns;
whichSlices= (min(kdx)+ round(0.5*(nslicesPerSlab))): round(nslicesPerSlab): max(kdx);
%whichSlices= whichSlices(2:end-1);
ns= length(whichSlices);
fovz= sliceLocations(whichSlices(end))- sliceLocations(whichSlices(1));

fovzHR= res* size(mask,3);

nskips=2;
maskMS= mask(1:nskips:end,1:nskips:end, whichSlices);
b1mapsMSn= b1maps(1:nskips:end,1:nskips:end,whichSlices,:);
b0mapMS= b0map(1:nskips:end,1:nskips:end,whichSlices);

%% manually fix the mask if needed
% if bDrawROI
%     maskMS1= maskMS;
%     figure;
%     for idx=1:size(maskMS1,3)
%         myimagesc(maskMS(:,:,idx));
%         imask= roipoly;
%         maskMS1(:,:,idx)= imask;
%     end
%     maskMS= maskMS1;
%     %figure, myMontagemn(maskMS,ns,1), colormap jet
% end
% %%%=================

foxktMS = 1e-3*[fov fovz]; % field of excitation in m
foxkt= 1e-3*[fov fovzHR]; 
poffset = [0 0 0]; % field of view shift (in mm) that has been applied to B1+ and B0 mapping. 

% water imaging
spect = [zeros(size(fatbnd)) ones(size(waterbnd))];
wts = [wsb*ones(size(fatbnd)) wpb*ones(size(waterbnd))];

% % fat imaging
% spect = [ones(size(fatbnd)) zeros(size(waterbnd))];
% wts = [wpb*ones(size(fatbnd)) wsb*ones(size(waterbnd))];

% use CP phase as a start point
nptsPerGrp= 3;
ngrps= 8;
x00= repmat(x0,[1 nptsPerGrp*ngrps]);
x00= x00.';

%%
spspktpointDesigner = SPSPkTpointPulseDesigner(b1mapsMSn,maskMS,foxktMS);
spspktpointDesigner.B0Map = b0mapMS;
spspktpointDesigner.DwellTime = dt;
spspktpointDesigner.MaxGradSlewRate = 100;%160;
spspktpointDesigner.Lambda = [1e-5 1e-4 1e-3 1e-2 1e-1];
spspktpointDesigner.NumOfPoints= nptsPerGrp; % current implementation only supports 3 kt point per group. 
spspktpointDesigner.NumOfGroups= ngrps;
spspktpointDesigner.x0= x00(:);

spspktpointDesigner.Frequencies= freqs;
spspktpointDesigner.Spectrum = spect;
spspktpointDesigner.Weights = wts;
%spspktpointDesigner.DesiredBandwidth=1400;
spspktpointDesigner.ReadOutOffset = poffset; % 
%spspktpointDesigner.ConvergenceTolerance= 1e-5;
%spspktpointDesigner.plotLCurve;

%
spspktpointDesigner.OptimalLambda = 1e-4;

[rf,grad]= spspktpointDesigner.design;
kp = spspktpointDesigner.Kspace;
rf= fa.* rf; 

gobj= gradPulse(grad,dt);
%gobj= gradPulse(switch_grad_polarity(grad,[1 1 -1]),dt);
figure, gobj.plot
figure, gobj.plotTraj
hold on
plot3(kp(1,:),kp(2,:),kp(3,:),'ro')
hold off

% rfobj= rfPulse(rf,dt);
% figure, rfobj.plot_amp
 
% save rfgrad rf grad

% high res
mxypatHR = run_bloch_sim (rf,grad,b1maps,[],foxkt,b0map,...
    0,[],dt,poffset); % a quick check on resulting |mxy| on resonance.

figure, orthosliceViewer(asind(abs(mxypatHR))), clim([0 15]), colormap jet
%figure, position_plots(asind(abs(mxypat)),[1 6],[],[],maskMS);
%save mxypatHR mxypatHR  


%% check frequency response
%freqs0=-500:25:500;% hz
freqs0= min(freqs):50:max(freqs);
%b0s= hz2tesla(repmat(freqs.',[1, size(rf,2)]));

clear mzpats3d mxypats3d
for idx=1:length(freqs0)
    
 imxypat = run_bloch_sim (rf,grad,b1mapsMSn,maskMS,foxkt,...
     b0mapMS+hz2tesla(freqs0(idx)),...
    0,[],dt,poffset);

    %mzpats3d(:,:,:,idx)= imzpatptx3d;
    mxypats3d(:,:,:,idx)= imxypat;
end

%
freqresp=zeros(1,size(mxypats3d,4));
for idx=1:size(mxypats3d,4)
    imxy= mxypats3d(:,:,:,idx);
    freqresp(idx)= mean(abs(imxy(maskMS)));
end

figure, plot(freqs0,(freqresp/max(freqresp)),'o-')
title('Frequency response')
xlabel('Frequency (Hz)')
ylabel('|Mxy| (a.u.)')

