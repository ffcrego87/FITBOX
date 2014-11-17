function msfcn_ff(block)
% Level-2 M file S-Function for unit delay demo.
%   Copyright 1990-2004 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $ 

  setup(block);
  
%endfunction

function setup(block)
  block.NumDialogPrms  = 40;
  
  %% Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 5;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  %block.InputPort(1).Dimensions        = -1;
  block.InputPort(1).DirectFeedthrough = true;
  
  block.InputPort(2).Dimensions        = -1;
  block.InputPort(2).DirectFeedthrough = true;
  
  block.OutputPort(1).Dimensions       = 1; %Fault 
  block.OutputPort(2).Dimensions       = block.DialogPrm(1).Data.nx; %upper bound
  block.OutputPort(3).Dimensions       = block.DialogPrm(1).Data.nx; %lower bound
  block.OutputPort(4).Dimensions       = 1; %mode
  block.OutputPort(5).Dimensions       = 1; %false_dtc
  
  %% Set block sample time to inherited
  Ts = block.DialogPrm(2).Data; %Sample time
  block.SampleTimes = [Ts 0]; %VERY IMPORTANT TO HAVE THE RIGHT SAMPLE TIME
  
  %% Register methods
  block.RegBlockMethod('SetInputPortSamplingMode',@SetInputPortSamplingMode);
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('SetInputPortDimensions',  @SetInputPortDimensions);  
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('InitializeConditions',    @InitConditions); 
  block.RegBlockMethod('Update',                  @Update);  
  %block.RegBlockMethod('SetOutputPortDimensions', @SetOutputPortDimensions);
  
%endfunction

function DoPostPropSetup(block)
  %% Setup Dwork
  block.NumDworks = 1;
  %controller variables
  block.Dwork(1).Name = 'k'; 
  block.Dwork(1).Dimensions      = 1;
  block.Dwork(1).DatatypeID      = 0;
  block.Dwork(1).Complexity      = 'Real';
  block.Dwork(1).UsedAsDiscState = true;
%endfunction

function InitConditions(block)

  %% Initialize Dwork
  block.Dwork(1).Data = 1;
%endfunction


function Update(block)
  k = block.Dwork(1).Data;
  block.Dwork(1).Data = k+1;
%endfunction

function Output(block)
    k = block.Dwork(1).Data;
    
    ff= block.DialogPrm(1).Data;
    Ts = block.DialogPrm(2).Data; %Sample time
    %Nominal
    A = block.DialogPrm(3).Data;
    B = block.DialogPrm(4).Data;
    L = block.DialogPrm(5).Data;
    C = block.DialogPrm(6).Data;
    D = block.DialogPrm(7).Data;
    N = block.DialogPrm(8).Data;
    wp= block.DialogPrm(9).Data;
    wm= block.DialogPrm(10).Data;
    %Fault3
    A3 = block.DialogPrm(11).Data;
    B3 = block.DialogPrm(12).Data;
    L3 = block.DialogPrm(13).Data;
    C3 = block.DialogPrm(14).Data;
    D3 = block.DialogPrm(15).Data;
    N3 = block.DialogPrm(16).Data;
    wp3 = block.DialogPrm(17).Data;
    wm3 = block.DialogPrm(18).Data;
    %Global
    Aglb = block.DialogPrm(19).Data;
    Bglb = block.DialogPrm(20).Data;
    Lglb = block.DialogPrm(21).Data;
    Cglb = block.DialogPrm(22).Data;
    Dglb = block.DialogPrm(23).Data;
    Nglb = block.DialogPrm(24).Data;
    %Fault2
    A2 = block.DialogPrm(25).Data;
    B2 = block.DialogPrm(26).Data;
    L2 = block.DialogPrm(27).Data;
    C2 = block.DialogPrm(28).Data;
    D2 = block.DialogPrm(29).Data;
    N2 = block.DialogPrm(30).Data;
    wp2 = block.DialogPrm(31).Data;
    wm2 = block.DialogPrm(32).Data;
    %Fault1
    A1 = block.DialogPrm(33).Data;
    B1 = block.DialogPrm(34).Data;
    L1 = block.DialogPrm(35).Data;
    C1 = block.DialogPrm(36).Data;
    D1 = block.DialogPrm(37).Data;
    N1 = block.DialogPrm(38).Data;
    wp1 = block.DialogPrm(39).Data;
    wm1 = block.DialogPrm(40).Data;
    %Faulty input data
    
    
    u = block.InputPort(1).Data;
    y = block.InputPort(2).Data;
    tic
    [id,active,false_dtc] = ff.update({'nominal'},A,B,L,C,D,N,u,y,wp,wm,...
        {'fault1'},A1,B1,L1,C1,D1,N1,u,y,wp1,wm1,...
        {'fault2'},A2,B2,L2,C2,D2,N2,u,y,wp2,wm2,...
        {'fault3'},(A3+A)/2,(B3+B)/2,(L3+L)/2,C3,D3,N3,u,y,wp3,wm3,...
        {'global'},Aglb,Bglb,Lglb,Cglb,Dglb,Nglb,[],[],[],[]);
    [x,d] = ff.getStateEstimate(id);
    disp(['k = ' num2str(k) ':' num2str(toc) ' s, active svos: ' num2str(ff.active)])
    
    idx = find(active,1);
    block.OutputPort(1).Data = idx; %fault number
    block.OutputPort(2).Data = x+d/2; %upper state bound
    block.OutputPort(3).Data = x-d/2; %lower state bound
    block.OutputPort(4).Data = ff.mode; %mode
    block.OutputPort(5).Data = false_dtc; %mode
%endfunction

function SetInputPortSamplingMode(block, idx, fd)
    block.InputPort(idx).SamplingMode = fd;
    block.InputPort(idx).SamplingMode = fd;
    block.OutputPort(1).SamplingMode = fd;
    block.OutputPort(2).SamplingMode = fd;
    block.OutputPort(3).SamplingMode = fd;
    block.OutputPort(4).SamplingMode = fd;
    block.OutputPort(5).SamplingMode = fd;
%endfunction

function SetInputPortDimensions(block, port, dimsInfo)
    block.InputPort(port).Dimensions = dimsInfo;
%endfunction
