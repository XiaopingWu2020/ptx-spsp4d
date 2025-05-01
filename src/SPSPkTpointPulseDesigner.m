%%% SPSPkTpointPulseDesigner.m --- 
%% 
%% Filename: SPSPkTpointPulseDesigner.m
%% Description: 
%% Author: Xiaoping Wu <xpwu@cmrr.umn.edu> 
%% Maintainer: 
%% Copyright (C) 2011 CMRR at UMN
%% Created: Fri Nov 18 11:05:14 2011 (CST)
%% Version: 
%% Last-Updated: Mon Dec 19 21:29:21 2011 (CST)
%%           By: Xiaoping Wu
%%     Update #: 249
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Commentary: 
%% 
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Change log:
%% 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
%%% Code:

classdef SPSPkTpointPulseDesigner < handle
    properties 
      B1Map = []
      Mask = []
      FieldOfExcit = [] 
      NumOfPoints = 3
      ConvergenceTolerance = 1e-5
      Lambda = [1e-1 1e0 1e1 1e2 1e3 1e4]
      OptimalLambda = 10
      B0Map = []
      MaxGradSlewRate = 160 
      MaxGradAmplitude = 50e-3
      DwellTime = 10e-6
      ReadOutOffset = [0 0 0] % mm
      Frequencies = 0
      Weights = 1
      Spectrum = 1
      NumOfGroups = 1
      DesiredBandwidth = 2*1050
      B1Constraint = inf
      ReshapeMethod = 'suppress'
      x0= 0;
    end % properties

    properties (SetAccess = protected)
    end % properties

    properties (SetAccess = private)
      IsNewGradient = true
      Gradient = {}
      Kspace = []
      BasicRF = {}
      IsNewBasicRF = true
    end % properties
 
    methods 
        function obj = SPSPkTpointPulseDesigner (b1map,mask,fox)
        % Constructor
        % Usage: obj = classname (a,b,c)
        % a:
        % b:
        % c:
        
        % Pre init.
        % anything not using output (obj)

          if nargin == 0
            disp(['-> The object is constructed. Please specify necessary ' ...
                  'properties...'])
            return
          end

        % compvalue = classname.staticMethod();

        
        % Object init.
        % call super class, if applicable, before accessing object
        % obj = obj@superClass(args{:});

        
        % Post init.
        % anything including accessing object
        % obj.classMethod();
        % obj.Property = compvalue;
          obj.B1Map = b1map;
          obj.Mask = mask;
          obj.FieldOfExcit = fox;
          
          disp('-> Object constructed. Please specify other needed properties...')     
        end
        
        function obj = set.ReshapeMethod (obj, rsmethod)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        if ~(strcmpi(rsmethod,'equalize')||strcmpi(rsmethod,'suppress'))
          disp('ReshapeMethod must be equalize or suppress')
          rsmethod = 'equalize';
        end
        
        obj.ReshapeMethod = rsmethod;

        end

        
        function obj = set.NumOfPoints (obj, npts)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
          if npts< 1
            npts = 1;
            disp('-> NumOfPoints must be >=1 and has been set to 1.')
          end
        
          if ~isequal(obj.NumOfPoints,npts)          
            obj.NumOfPoints = npts;          
            obj.IsNewGradient = true;
            obj.IsNewBasicRF = true;
          end        
        
        end
        
        function obj = set.NumOfGroups (obj, ngrps)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        if ngrps< 1
          ngrps = 1;
          disp('-> NumOfGroups must be >=1 and has been set to 1.')
        end

        if ~isequal(obj.NumOfGroups,ngrps)          
          obj.NumOfGroups = ngrps;          
          obj.IsNewGradient = true;
          obj.IsNewBasicRF = true;
        end        
          
        end


        function obj = set.B0Map (obj, b0map)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        if isempty(b0map)
          b0map = zeros(size(obj.Mask));
        end  
        
        if ~isequal(obj.B0Map, b0map)
          obj.B0Map = b0map;
          obj.IsNewBasicRF = true;
        end

        end
        
        function obj = set.Frequencies (obj, freqs)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        
        if ~isequal(obj.Frequencies, freqs)
          obj.Frequencies = freqs;
          obj.Spectrum = ones(size(freqs));
          obj.Weights = ones(size(freqs));
          obj.IsNewBasicRF = true;
        end
          
        end
        
        function obj = set.Spectrum (obj, spect)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        
        if length(spect)~= length(obj.Frequencies)
          error('-> Spectrum size mismatches with Frequencies size...')
        end
        
        if ~isequal(obj.Spectrum, spect)
          obj.Spectrum = spect;
          obj.IsNewBasicRF = true;
        end
        
        end
          
        function obj = set.OptimalLambda (obj, optLam)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        
        if obj.OptimalLambda ~= optLam
          obj.OptimalLambda = optLam;
          obj.IsNewBasicRF = true;
        end
        
        end
        
        function obj = set.Weights (obj, wts)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        
        if length(wts)~= length(obj.Frequencies)
          error('-> Spectrum size mismatches with Frequencies size...')
        end
        
        if ~isequal(obj.Weights, wts)
          obj.Weights = wts;
          obj.IsNewBasicRF = true;
        end
          
        end


        function stat = plotLCurve (obj)
        % Purpose:
        % Usage: stat = obj.plotLCurve ()
        % 
        % 
        % 
        lambda = obj.Lambda;
        obj.designKspace;
        [~,rho,eta] = obj.calcWeights([],lambda,2);
        figure, plot_lc(rho,eta,'o',1,lambda)
        
        stat = 1;
        end

        function [myrf,mygrad] = design (obj)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:

        if obj.IsNewGradient
          obj.designKspace;
          obj.designGrad;
          obj.IsNewGradient = false;
        end

        if obj.IsNewBasicRF
          subrfLen = obj.calcSubRfLengthBasic;
          phasetrack = obj.calcPhasetrack(subrfLen);
          rf = obj.calcRF(phasetrack,subrfLen);
        end          
        
        % equalize RF amp 
        if obj.hasHighRfPeak

          switch obj.ReshapeMethod
           case 'equalize'
            subrfLen = obj.calcSubRfLengthEqu;
           case 'suppress'
            subrfLen = obj.calcSubRfLengthSup;            
           otherwise
          end

          phasetrack = obj.calcPhasetrack(subrfLen);
          rf = obj.calcRF(phasetrack,subrfLen);
        end
                
        [myrf,mygrad] = obj.assemblePulse(rf);
        
        disp('-> kT points SPSP pulse design DONE...')
          
        end

        
        function subrfLen = calcSubRfLengthSup (obj)
          % only suppress peak rf amp
          rf = obj.BasicRF;
          subrfLen = zeros(size(rf));
          for ind= 1:length(subrfLen),
            irf = rf{ind};
            ampmax = max(abs(irf(:)));
            
            if ampmax > obj.B1Constraint
              subrfLen(ind) = ceil(size(irf,2).* ampmax./ obj.B1Constraint);
            else
              subrfLen(ind) = size(irf,2);
            end
           
          end
          
        end
     
        function subrfLen = calcSubRfLengthEqu (obj)
          % equalize or VERSE rf amp
        
          rf = obj.BasicRF;
          subrfLen = zeros(size(rf));
          for ind= 1:length(subrfLen),
            irf = rf{ind};
            subrfLen(ind) = ceil(size(irf,2).* max(abs(irf(:)))./ obj.B1Constraint);
          end
        end
        
        function isHigh = hasHighRfPeak (obj)
          rf = obj.BasicRF;
          ampmax = 0;
          for ind= 1: length(rf),
            ampmax = max([ampmax; max(abs(rf{ind}(:)))]);            
          end
          
          if ampmax> obj.B1Constraint
            isHigh = true;else
            isHigh = false;
          end
          
        end
        
        function rf = calcRF (obj, phasetrack, subrfLen)
        
        [wt1,rho1] = obj.calcWeights(phasetrack,obj.OptimalLambda,1);
        wt11 = reshape(wt1,[],size(obj.B1Map,4));
        wt11 = wt11.';
  
        rf = cell(size(subrfLen));
        for ind=1:length(subrfLen),
          nt = subrfLen(ind);
          rf{ind} = 1/nt.* wt11(:,ind) * ones(1,nt);
        end

        if obj.IsNewBasicRF
          obj.BasicRF = rf;
          obj.IsNewBasicRF = false;
        end
        
        end
                
        function [myrf, mygrad] = assemblePulse (obj, rf)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        grad = obj.Gradient;
        mygrad = [];
        myrf = [];
        nchs = size(rf{1},1);
        npts = length(rf);
        for ind=1:npts,
          mygrad = [mygrad zeros(3,size(rf{ind},2)),grad{ind}];
          myrf = [myrf rf{ind} complex(zeros(nchs,size(grad{ind},2)))];
        end
        
        mygrad = [mygrad zeros(3,1)];
        myrf(:,size(mygrad,2)) = 0;
        
        end

      function designKspace (obj)
      % Purpose:
      % Usage: output = obj.funcname (a,b,c)
      % a:
      % b:
      % c:
      
%       % new implementation which did not seem to work well since it
%       % resulted in rf energy deposited solely at the kspace center. 
%       % 
%       ngrps = obj.NumOfGroups;
%       fox= obj.FieldOfExcit;
%       myfox= mean(fox);
%       if fox(3)/mean(fox(1:2))< 0.1 % inplane sampling only
%           mytheta0 = 0:180/ngrps:180;
%           mytheta= mytheta0;
%           kp0=[0 0 0;1 0 0;-1 0 0].';
%           
%           kp= [];
%           for ind=1:length(mytheta),
%               ikp = angle2dcm(deg2rad(mytheta(ind)),0,0,'ZYX') * kp0;
%               kp = [kp ikp];
%           end
%         
%       else
%           mytheta0 = 0:180/ngrps:180;
%           %
%           mytheta= mytheta0(1:2:end); % put half in the inplane
%           kp0=[0 0 0;1 0 0;-1 0 0].';
%           
%           kp= [];
%           for ind=1:length(mytheta),
%               ikp = angle2dcm(deg2rad(mytheta(ind)),0,0,'ZYX') * kp0;
%               kp = [kp ikp];
%           end
%           %
%           myphi= 45*fox(3)/mean(fox(1:2)); % deg
%           mytheta1= mytheta0(2:2:end);
%           
%           % plus phi
%           mytheta11= mytheta1(1:2:end);
%           kp0=[0 0 0;...
%               cosd(myphi) 0 sind(myphi);...
%               cosd(myphi+180) 0 sind(myphi+180)].';
%          
%           for ind=1:length(mytheta11),
%               ikp = angle2dcm(deg2rad(mytheta11(ind)),0,0,'ZYX') * kp0;
%               kp = [kp ikp];
%           end
%           % minus phi
%           mytheta12= mytheta1(2:2:end);
%           kp0=[0 0 0;...
%               cosd(-myphi) 0 sind(-myphi);...
%               cosd(-myphi+180) 0 sind(-myphi+180)].';
%           
%           for ind=1:length(mytheta12),
%               ikp = angle2dcm(deg2rad(mytheta12(ind)),0,0,'ZYX') * kp0;
%               kp = [kp ikp];
%           end
%         
%           
%       end
% 
%       obj.Kspace = 2*pi./myfox.* kp;
      
      % old implementation which did not present a uniform sampling of
      % kspace.
      ngrps = obj.NumOfGroups;
      myfox = mean(obj.FieldOfExcit);
      mytheta = 0:180/ngrps:180;
      mytheta = mytheta(1:end-1);
      myphi = mytheta;

      switch obj.NumOfPoints
       case 3
        kp0 = [0 0 0;
              -1 0 0;
              1 0 0].';
          otherwise
           error('-> NumOfPoints not supported...')
      end
      
      kp= []; 
      dk= 0.25*pi./myfox;
      for ind=1:length(mytheta),
        ikp = angle2dcm(deg2rad(mytheta(ind)),deg2rad(myphi(ind)),0,'ZYX') * kp0;
        kp = [kp ind*dk*ikp];%
        %kp = [kp ceil(ind/2)*dk*ikp];
      end
           
      obj.Kspace = kp; %2*pi./myfox.* kp;
      
      end


       function phasetrack = calcPhasetrack (obj,subrfLen)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
          grad = obj.Gradient;
          phasetrack= subrfLen;
          for ind= 1: length(grad)-1,
            phasetrack(ind+1) = phasetrack(ind)+ size(grad{ind},2)+ subrfLen(ind+1);
          end
        end
  
        function subrfLen = calcSubRfLengthBasic (obj)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:
        grad = obj.Gradient;
        npts = length(grad);
        subrfLen = zeros(1,npts);

        bw= obj.DesiredBandwidth;
        dt = obj.DwellTime;
        dt0 = 1/bw;
        ntimepts = round(dt0/dt);
        
        npts0= obj.NumOfPoints;
        
        % calc sub rf lengths
        for ind=1:npts0:npts,
          g0=[];
          for jdx = 1:npts0
            g0= [g0 grad{ind+jdx-1}];
          end
          
          nt = round((ntimepts - size(g0,2))/npts0);          
          if nt>0
            subrfLen(ind:ind+npts0-1) = nt; else
            subrfLen(ind:ind+npts0-1) = 1;
          end
        end
        
        end


        function designGrad (obj)
        % Purpose:
        % Usage: output = obj.funcname (a,b,c)
        % a:
        % b:
        % c:

        if isempty(obj.Kspace)
          obj.designKspace;
        end
        
        maxamp = obj.MaxGradAmplitude;
        maxsr = obj.MaxGradSlewRate;
        dt = obj.DwellTime;
        
        kp = [obj.Kspace zeros(3,1)]; % back to origin
        dkval = diff(kp,1,2);
        grad = {};
        nblips = size(dkval,2);
        for ind= 1: nblips     
          dkmax = max(abs(dkval(:,ind)));
          %blipmax = design_toptgrad1D(0,0,dkmax,maxamp,maxsr,dt);
          blipmax = design_grad_trapz(dkmax,maxamp,maxsr,dt);
          blips = [blipmax;blipmax;blipmax];
          sf = diag(dkval(:,ind)./ dkmax);
    
          grad{ind} = sf*blips;  
        end
        obj.Gradient = grad;
          
        end

      function [wt,rho,eta] = calcWeights (obj,phasetrack,lambda,itlv)
      % Purpose:
      % Usage: output = obj.funcname (a,b,c)
      % a:
      % b:
      % c:

      if isempty(obj.Kspace)
        obj.designKspace;
      end
      
        kp = obj.Kspace;
        dt = obj.DwellTime;
        b1map = obj.B1Map;
        mask = obj.Mask;
        fox = obj.FieldOfExcit;
        b0map = obj.B0Map;
        if isempty(phasetrack)
          b0map = [];
        end
        

        %phas= -pi: 2*pi/16:pi

        freqs = obj.Frequencies(1:itlv:end);
        wts = obj.Weights(1:itlv:end);
        spect = obj.Spectrum(1:itlv:end);
        poffset = obj.ReadOutOffset;

        m = construct_targvect_spsp(obj.Mask,obj.Mask,spect);

        sysmat = construct_sysmat_spspkT(kp,b1map,mask,fox,b0map, phasetrack, ...
                                         freqs,wts,dt,poffset);
        
        tol = obj.ConvergenceTolerance;
        
        [wt,rho,eta]= solve_mlstr(sysmat,m,lambda,tol,obj.x0);

      end

    end % methods
end % classdef



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SPSPkTpointPulseDesigner.m ends here
