classdef fim < handle
%
% fim - Fault Identification Module
%
% Requirements: MATLAB 2011a 
%               Parallel Computing Toolbox
%
% Description: This class manages the set valued estimates of a bank of
% SVOs in order to identify fault events.
%
% Properties:
%   SVObank - (Protected) array of SVOs. 
%             Use method 'addSVO' to increase the number of SVOs. 
%             Must have a nominal SVO and a global SVO.
%   mode    - (Protected)
%             Current operating mode of the fault identification module.
%   nx      - (Immutable) number of the state variables.
%   nu      - (Immutable) number of inputs.
%   ny      - (Immutable) number of outputs.
%   nw      - (Immutable) number of disturbances.
%   x0      - (Immutable) estimate of the initial state
%   unc     - (Immutable) uncertainty in x0
%   Treset  - (Immutable) timer limit for fault identification
%   active  - (Protected) set of active SVOs
%   ids     - (Protected) set of SVO identification strings
%   classes - (Protected) array that stores each SVO class
% Methods
%   addSVO(id,class,n,options{approximation,solver,angTol,numTol})
%   [fault_no,activeSVOs,false_dtc] = update(A,B,L,C,D,N,u,y,wp,wm)
%
%     Copyright (C) 2013  Pedro Casau
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%     
    properties (SetAccess = immutable)
        nx;
        Treset;
    end
    properties (SetAccess = protected)
        SVObank;
        mode;
        active;
        ids;
        classes;
    end
    properties (Hidden = true)
        update_flag;
        nSVO;
        nml;
        glb;
        timer;
        numWorkers;
    end
    methods
        function obj = fim(Treset,nx)
            %fim constructor
            if nargin >0
                if nargin ~= 2
                    ME = MException('fim:constructor:invalidInput',...
                        'Constructor of type FIM must have 2 input arguments.');
                    throw(ME);
                elseif ~isscalar(Treset) || Treset ~= uint64(Treset) || Treset<=0
                        ME = MException('fim:constructor:invalidInput',...
                            'Treset must be a scalar greater than 0.');
                        throw(ME);
                elseif ~isscalar(nx) || nx ~= uint64(nx) || nx <= 0
                        ME = MException('fim:constructor:invalidInput',...
                            'nx must be a scalar greater than 0.');
                        throw(ME);
                end
                obj.nx      = nx;
                obj.SVObank = [];
                obj.ids     = cell(0);
                obj.classes = [];
                obj.active  = [];
                obj.update_flag = false;
                obj.mode    = 1; %nominal mode is the default mode
                obj.nSVO    = 0;
                obj.glb     = 0;
                obj.nml     = 0;
                obj.Treset  = Treset;
                obj.numWorkers = 0;
            end
        end
        
        function addSVO(obj,id,class,svos)
            if obj.update_flag
                ME = MException('fim:addSVO:invalidAction',...
                    'Cannot add SVO after a call to update has been made.');
                throw(ME);
            elseif nargin < 3
                ME = MException('fim:addSVO:notEnoughInputs',...
                    'Not enough input arguments.');
                throw(ME);
            elseif all([1 2] ~= class)
                ME = MException('fim:addSVO:invalidInput',...
                    'class must be either 1 or 2');
                throw(ME);
            elseif ~ischar(id) || any(strcmpi(obj.ids,id))
                ME = MException('fim:addSVO:invalidInput',...
                    'id must be a unique string');
                throw(ME);
            elseif ~isa(svos,'SVO')
                ME = MException('fim:addSVO:invalidInput',...
                    'Last input must belong to the class "SVO"');
                throw(ME);
            elseif svos.nx ~= obj.nx
                ME = MException('fim:addSVO:invalidInput',...
                    'The given SVO must have a number of states equal to %d',obj.nx);
                throw(ME);
            end
            if strcmpi(id,'nominal')
                if class ~= 1
                    ME = MException('fim:addSVO:invalidInput',...
                        'Nominal SVO must be class 1.');
                    throw(ME);
                end
                obj.nml = obj.nSVO+1;
            elseif strcmpi(id,'global')
                if class ~= 1
                    ME = MException('fim:addSVO:invalidInput',...
                        'Global SVO must be class 1.');
                    throw(ME);
                end
                obj.glb = obj.nSVO+1;
                %options.approximation = 'hypercube'; %this property is used 
                    %in fim:update during the nominal mode operation
            end
            obj.nSVO = obj.nSVO + 1;
            obj.SVObank = [obj.SVObank;
                svos];
            obj.ids = [obj.ids;id];
            obj.classes = [obj.classes;class];
            if class == 1
                %Class 1 SVOs are active by default
                obj.active = [obj.active;1];
            else
                obj.active = [obj.active;0];
            end
        end
        
        function [fault_id,activeSVOs,false_dtc] = update(obj,varargin)
            if ~obj.nml
                ME = MException('fim:update:missingSVO',...
                    'Nominal SVO is not defined.');
                throw(ME);
            elseif ~obj.glb
                ME = MException('fim:update:missingSVO',...
                    'Global SVO is not defined.');
                throw(ME);
            end
            idx = ones(obj.nSVO,1);
            for I = 1:11:nargin-1
                if ~iscellstr(varargin{I})
                    ME = MException('fim:update:wrongInput',...
                        'Invalid input arguments. The user must supply sets of arguments ({ids},A,B,L,C,D,N,u,y,wp,wm).');
                    throw(ME);
                end
                for J = 1:numel(varargin{I})
                    index = find(strcmp(varargin{I}(J),obj.ids),1);
                    if isempty(index)
                        ME = MException('fim:update:wrongInput',...
                            'Invalid input arguments. The SVO id %s does not exist.',varargin{I}(J));
                        throw(ME);
                    end
                    idx(index) = 0;
                    A = varargin{I+1};
                    B = varargin{I+2};
                    L = varargin{I+3};
                    C = varargin{I+4};
                    D = varargin{I+5};
                    N = varargin{I+6};
                    u = varargin{I+7};
                    y = varargin{I+8};
                    wp= varargin{I+9};
                    wm= varargin{I+10};
                    save(['SVOdata' num2str(index)],'A','B','L','C','D','N','u','y','wp','wm')
                end
            end
            if any(idx)
                ME = MException('fim:update:missingData',...
                    'The user did not specify update data for some SVO.');
                throw(ME);
            end
                
                
            if ~obj.update_flag
                %prevents further access to addSVO
                obj.update_flag = true;
                %%%NEWER CODE%%%
                %schd = parcluster(parallel.defaultClusterProfile);
                %obj.numWorkers = schd.NumWorkers;
                %%if available, 1 worker per SVO is selected
                %matlabpool(parallel.defaultClusterProfile,min(obj.numWorkers,obj.nSVO));
                %%%%%%%%%%%%%%%%
                %%%CODE FOR MATLAB 2011a Compatibility%%%
                schd = findResource('scheduler', 'configuration', 'local');
                obj.numWorkers = schd.ClusterSize;
                %if available, 1 worker per SVO is selected
                matlabpool('local',min(obj.numWorkers,obj.nSVO));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            %%%%%%%%%%%%%% PARALLEL COMPUTING ROUTINE %%%%%%%%%%%%%%
            %save data
            for I = 1:obj.nSVO
                svo = obj.SVObank(I);
                a   = obj.active(I);
                save(['svo' num2str(I)],'svo');
                save(['act' num2str(I)],'a');
            end
            %run parallel for loop
            for I = 1:obj.nSVO
                aux = load(['svo' num2str(I)]);
                svos(I) = SVO();
                svos(I) = aux.svo;
                aux = load(['act' num2str(I)]);
                act(I) = aux.a;
                aux = load(['SVOdata' num2str(I)]);
                A = aux.A;
                B = aux.B;
                L = aux.L;
                C = aux.C;
                D = aux.D;
                N = aux.N;
                u = aux.u;
                y = aux.y;
                wp= aux.wp;
                wm= aux.wm;                
                if act(I)
                    act(I) = ~svos(I).update(A,B,L,C,D,N,u,y,wp,wm);
                end
            end
            %load data from distributed workers
            obj.active = act;
            obj.SVObank = svos;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ~obj.active(obj.glb)
                ME = MException('fim:update:globalSVO',...
                    'The set valued estimate of the global SVO is empty');
                throw(ME);
            end
            if obj.mode == 1 && ~obj.active(obj.nml) %Nominal mode
                [x,du] = obj.getStateEstimate('global');
                %%%%%%%%%%%%%% PARALLEL COMPUTING ROUTINE %%%%%%%%%%%%%% 
                %save data
                for I = 1:obj.nSVO
                    svo = obj.SVObank(I);
                    a   = obj.active(I);
                    clss= obj.classes(I);
                    save(['svo' num2str(I)],'svo');
                    save(['act' num2str(I)],'a');
                    save(['clss' num2str(I)],'clss');
                end
                %run parfor
                for I = 1:obj.nSVO
                    aux = load(['svo' num2str(I)]);
                    svos(I) = SVO();
                    svos(I) = aux.svo;
                    aux = load(['act' num2str(I)]);
                    act(I) = aux.a;
                    aux = load(['clss' num2str(I)]);
                    clss(I) = aux.clss;
                    if ~act(I) && clss(I) == 2
                        aux = load(['SVOdata' num2str(I)]);
                        A = aux.A;
                        B = aux.B;
                        L = aux.L;
                        C = aux.C;
                        D = aux.D;
                        N = aux.N;
                        u = aux.u;
                        y = aux.y;
                        wp= aux.wp;
                        wm= aux.wm;   
                        svos(I).reset(x,du);
                        act(I) = ~svos(I).update(A,B,L,C,D,N,u,y,wp,wm);
                    end
                end
                obj.active = act;
                obj.SVObank = svos;
                obj.mode = 2;
                obj.timer= 0;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            
            if sum(obj.active) == 1 %only global SVO is not faulty
                false_dtc = 1;
                reset = 1;
                index = obj.glb;
            else
                index = find([obj.active(1:obj.glb-1) 0 obj.active(obj.glb+1:end)],1);
                reset = 0;
                false_dtc = 0;
            end
            activeSVOs = obj.active;
            fault_id = obj.ids{index};
            if obj.mode == 2 && ~false_dtc
                if sum(obj.active) == 2
                    obj.mode = 3; %fault isolated
                elseif obj.timer >= obj.Treset
                    reset = 1;
                else
                    obj.timer = obj.timer +1;
                end
            end
            
            if obj.mode == 3
                if sum(obj.active) == 1 && obj.classes(index) == 1
                    false_dtc = 1;
                    reset = 1;
                elseif sum(obj.active) == 1 || ...
                        (obj.timer >= obj.Treset && obj.classes(index) == 1)
                    reset = 1;
                else
                    obj.timer = obj.timer + 1;
                end
            end
            
            if reset
                obj.mode = 1; %return to nominal
                [x,du] = obj.getStateEstimate(fault_id);
                %%%%%%%%%%%%%% PARALLEL COMPUTING ROUTINE %%%%%%%%%%%%%%
                %save data
                g = obj.glb;
                for I = 1:obj.nSVO
                    svo = obj.SVObank(I);
                    a   = obj.active(I);
                    clss= obj.classes(I);
                    save(['svo' num2str(I)],'svo');
                    save(['act' num2str(I)],'a');
                    save(['clss' num2str(I)],'clss');
                end
                %run parfor
                for I = 1:obj.nSVO
                        aux = load(['svo' num2str(I)]);
                        svos(I) = SVO();
                        svos(I) = aux.svo;
                        aux = load(['act' num2str(I)]);
                        act(I) = aux.a;
                        aux = load(['clss' num2str(I)]);
                        clss(I) = aux.clss;
                        if clss(I) == 1
                            %resets class 1 SVOs
                            aux = load(['SVOdata' num2str(I)]);
                            A = aux.A;
                            B = aux.B;
                            L = aux.L;
                            C = aux.C;
                            D = aux.D;
                            N = aux.N;
                            u = aux.u;
                            y = aux.y;
                            wp= aux.wp;
                            wm= aux.wm;   
                            svos(I).reset(x,du);
                            act(I) = ~svos(I).update(A,B,L,C,D,N,u,y,wp,wm);
                        else
                            act(I) = 0;
                        end
                end
                obj.active = act;
                obj.SVObank = svos;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
        end
        
        function delete(obj)
            %class destructor
            if obj.numWorkers ~= 0
                matlabpool close
            end
        end
    %end
    %methods %(Access = private)
        function [x,u] = getStateEstimate(obj,id)
            %SVO must have the property 'approximation' set to 'hypercube'
            index = strcmpi(obj.ids,id);
            [x,u] = obj.SVObank(index).getStateEstimate();       
        end
    end
end