classdef SVO < handle
%
% File: Definition of the SVO class
% Theory:
%     The Set-Valued Observers (SVO) are mathematical objects that provide
%     set valued estimates of the state of a Linear Time Varying (LTV) system
%     of the form
%     
%     x[k+1] = A[k]x[k]+Bu[k]+Lw[k]
%     y[k] = Cx[k]+Du[k]+Nw[k]
%     
%     where x[k] denotes the state of the dynamic system, u[k] is the input
%     of the system, y[k] denotes the output of the system and w[k] comprises
%     both the disturbances and noise acting on the system. Valid set valued
%     estimates are obtained under the assumptions:
%       1) the distubances and sensor noise are bounded or, more precisely,
%     they belong to a given hypercube B,
%       2) the initial state of the system belongs to a given compact set
%     X(0).
%     Moreover, the set-valued estimate does not ground unbounded if the
%     system is observable (c.f.~\cite{Shamma1999}). If these assumptions are
%     not met, then the set-valued estimate might not be valid.
%     
%     The SVO algorithm is described in~\cite{Shamma1999} and is as follows.
%     Let W(y,k) and X(y,k) denote the set of all possible disturbances
%     and state values, respectively, at the time step k given a particular
%     sequence of outputs y(0), y(1), ..., y(k). Let \tilde{X}(y) denote
%     the set of all possible state values given a particular y and w in B.
%     
%     Initialization: X(y,0) = \tilde{X}(y(0))\cap X(0) =Set(M(0),m(0))
%     
%     Propagation: Realize the Fourier-Motzkin elimination of the polytope
%     P [x(k),w(k),x(k-1),w(k-1)]' <= p, with
%     
%     P = [    I,    0, -A(k-1), -L(k-1);
%             -I,    0,  A(k-1),  L(k-1);
%           C(k), N(k),       0,       0;
%          -C(k),-N(k),       0,       0;
%              0,    I,       0,       0;
%              0,   -I,       0,       0;
%              0,    0,       0,       I;
%              0,    0,       0,      -I;
%              0,    0,  M(k-1),       0];
%     
%     p = [B(k-1)u(k-1);
%          -B(k-1)u(k-1);
%          y(k);
%          -y(k);
%          w+(k);
%          -w-(k);
%          w+(k-1);
%          -w-(k-1)];
%     
%     where w+[k] and w-[k] are the upper and lower bounds defining the
%     hypercube B and (M(k),m(k)) = projection(P,p,idx), where
%     projection() is a function which performs the fourier motzkin
%     elimination of all but the variables indexed by idx={1,...,n_x}.
% Properties:
%     n     - (Immutable) horizon (default: n=1)
%     (M,m) - Current polytope
%     nx    - (Immutable)  number of states
%     nu    - (Immutable)  number of inputs
%     ny    - (Immutable)  number of outputs
%     nw    - (Immutable)  number of disturbances
%
% Methods:
%     obj = SVO() - constructor
%     update(A(k),B(k),L(k),C(k),D(k),N(k),u(k),y(k),w+(k),w-(k)) - updates
%     (M,m)
%     projection() - projection of Set(P,p) onto the dimensions
%     indexed by idx={1,...,n_x}.
%     idx={1,...,n_x}
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
        n; %Horizon can only be defined in constructor (default = 1)
        nx;%number of states in the dynamic system 
        nu;%number of control inputs
        ny;%number of outputs
        nw;%number of disturbances
        approximation; %options - none, hypercube (default = none)
        solver; %options: MATLAB, CPLEX (default = matlab)
        angTol; %angle tolerance (default = 1e-2)
        numTol; %absolute numeric tolerance (default = 1e-6)
        accelerator; %faster solver for n >= 5 (default = 'off')
    end
    properties (SetAccess = protected)
        %user may retrieve the polytope (M,m) but must use the method
        %projection() to compute it
        vol;
        M;
        m;
        Ad;
        Bd;
        Ld;
        Cd;
        Dd;
        Nd;
    end
    properties (Hidden = true)
        k;  %integer in the range 0-n used to keep track of the current
        %iteration
        A;B;L;u;  %memory variables that store the previous iteration
        %data
        flag_Delta;
        proxy_M;
        proxy_m;
        nAd;
        nBd;
        nLd;
        nCd;
        nDd;
        nNd;
    end
    
        methods %(Access = private)
        function [x,fval,exitflag] = LPsolve(obj,f,A,b)
            if strcmpi('cplex',obj.solver)
                try
                    [fval,x,exitflag] = lps(-f',A,b);
                catch ME
                    ME1 = MException('SVO:isempty:solverNotAvailable',...
                        'CPLEX solver failed to execute');
                    ME1 = addCause(ME1,ME);
                    throw(ME1);
                end
            else
                try
                    options = optimset('Display','off','MaxIter',1000);
                    [x,fval,exitflag] = linprog(f,A,b,[],[],[],[],[],options);
                catch ME
                    ME1 = MException('SVO:isempty:solverNotAvailable',...
                        'LINPROG solver failed to execute');
                    ME1 = addCause(ME1,ME);
                    throw(ME1);
                end    
            end
        end
        function [M,empty] = boundingBox(obj,P,p,idx)
            %Project Mx<=K onto dimensions idx
            j = size(P,2);
            l = numel(idx);
            aux = eye(j);
            upr = zeros(l,1);
            low = zeros(l,1);
            for I=1:l
                [x_low,~,status_low] = obj.LPsolve(aux(idx(I),:),P,p);
                [x_upr,~,status_upr] = obj.LPsolve(-aux(idx(I),:),P,p);
                upr(I) = x_upr(idx(I));
                low(I) = x_low(idx(I));
                if all(status_upr ~= [1 5]) || all(status_low ~= [1 5])
                    warning(['The set will be considered empty. Exitflags: ' num2str([status_upr status_low])]); 
                    empty = true;
                    M = [];
                    return
                end
            end
            empty = false;
            M = [eye(l),upr;
                -eye(l),-low];
            obj.vol = prod(upr-low);
        end
        function [H,h,empty] = cvxHull(obj,varargin)
            if nargin < 2
                ME = MException('SVO:cvxHull:notEnoughInputs',...
                    'Not enough input arguments');
                throw(ME)
            end
            for I = 1:numel(varargin)-1
                if size(varargin{I},2) ~= size(varargin{I+1},2)
                    ME = MException('SVO:cvxHull:inconsistentDimensions',...
                        'Each polytope must have the same number of columns');
                    throw(ME);
                end
            end
            nvar = size(varargin{1},2)-1;
            if strcmpi(obj.approximation,'none')
                ME = MException('SVO:cvxHull:underConstruction',...
                    'Exact convex hull implementation missing');
                throw(ME);
            else
                ub = -Inf*ones(nvar,1);
                lb = Inf*ones(nvar,1);
                for I = 1:numel(varargin)
                    P = varargin{I};
                    P = obj.boundingBox(P(:,1:end-1),P(:,end),1:nvar);
                    ub= max([ub,P(1:nvar,end)],[],2);
                    lb= min([lb,-P(nvar+1:2*nvar,end)],[],2);
                end
                if any(ub < lb)
                    empty = true;
                    H = [];
                    h = [];
                    warning('SVO:cvxHull:emptyPolytope','The polytope is empty');
                else
                    H = [eye(nvar);-eye(nvar)];
                    h = [ub;-lb];
                end
            end
        end
    end
    
    methods
        function obj = SVO(nx,nu,ny,nw,x0,uncertainty,n,options)
            %SVO object constructor
            if nargin > 0 %no input case generates an empty class instance
                if nargin < 6
                    ME = MException('SVO:constructor:notEnoughInputs',...
                        'Not enough input arguments for constructor of type SVO.');
                    throw(ME);
                elseif nargin > 8
                    ME = MException('SVO:constructor:tooManyInputs',...
                        'Too many input arguments for constructor of type SVO.');
                    throw(ME);
                elseif nargin >= 7
                    if ~isscalar(n) || n ~= uint8(n) || n<=0
                        ME = MException('SVO:constructor:invalidDimensions',...
                            'The SVO horizon n must be a natural number.');
                        throw(ME);
                    end
                    obj.n = n;   
                else
                    obj.n = 1;
                end
                if nargin == 8 && isstruct(options)
                    if ~isfield(options,'approximation') ||...
                            ~xor(strcmpi(options.approximation,'none'),...
                            strcmpi(options.approximation,'hypercube'))
                        ME = MException('SVO:constructor:invalidOptions',...
                            'Options must contain the field approximation, which is either hypercube or none');
                        throw(ME);
                    elseif ~isfield(options,'solver') ||...
                            ~xor(strcmpi(options.solver,'matlab'),...
                            strcmpi(options.solver,'cplex'))
                        ME = MException('SVO:constructor:invalidOptions',...
                            'Options must contain the field solver, which is either matlab or cplex');
                        throw(ME);
                    elseif ~isfield(options,'accelerator') ||...
                        ~xor(strcmpi(options.accelerator,'on'),...
                            strcmpi(options.accelerator,'off'))
                        ME = MException('SVO:constructor:invalidOptions',...
                            'Options must contain the field accelerator ("on" or "off").');
                        throw(ME);
                    end
                    obj.approximation = options.approximation;
                    obj.solver = options.solver;
                    obj.accelerator = options.accelerator;
                    if ~isfield(options,'angTol') ||...
                            ~isscalar(options.angTol) ||...
                            options.angTol <= 0
                        ME = MException('SVO:constructor:invalidOptions',...
                            'Options must contain the field angTol, which must be a scalar greater than 0');
                        throw(ME);
                    elseif ~isfield(options,'numTol') ||...
                            ~isscalar(options.numTol) ||...
                            options.numTol <= 0
                        ME = MException('SVO:constructor:invalidOptions',...
                            'Options must contain the field numTol, which must be a scalar greater than 0');
                        throw(ME);
                    end
                    %angular tolerance saturates at 10???
                    obj.angTol = min(options.angTol,10);
                    %absolute numeric tolerance saturates at 1
                    obj.numTol = min(options.numTol,1);
                else
                    obj.approximation = 'hypercube';
                    obj.angTol = 1e-2;
                    obj.numTol = 1e-6;
                    obj.solver = 'matlab';
                    obj.accelerator = 'off';
                end
                if ~isvector(x0)
                    ME = MException('SVO:constructor:invalidDimensions',...
                        'Initial conditions x0 must be a vector');
                    throw(ME);
                elseif ~isvector(uncertainty) || numel(x0) ~= numel(uncertainty)
                    ME = MException('SVO:constructor:inconsistentDimensions',...
                        'The uncertainty vector does not have the same dimensions as the initial conditions x0');
                    throw(ME);
                elseif nx ~= numel(x0)
                    ME = MException('SVO:constructor:inconsistentDimensions',...
                        'The number of state variables nx does not match the number of elements in the initial state x0.');
                    throw(ME);
                elseif ~isscalar(nx) || nx ~= uint8(nx) || nx<=0
                        ME = MException('fim:constructor:invalidInput',...
                            'nx must be a scalar greater than 0.');
                        throw(ME);
                elseif ~isscalar(nu) || nu ~= uint8(nu) || nu<0
                        ME = MException('fim:constructor:invalidInput',...
                            'nu must be a scalar greater or equal to 0.');
                        throw(ME);
                elseif ~isscalar(ny) || ny ~= uint8(ny) || ny<0
                        ME = MException('fim:constructor:invalidInput',...
                            'ny must be a scalar greater or equal to 0.');
                        throw(ME);
                elseif ~isscalar(nw) || nw ~= uint8(nw) || nw<0
                        ME = MException('fim:constructor:invalidInput',...
                            'nw must be a scalar greater or equal to 0.');
                        throw(ME);
                end
                
                x0 = reshape(x0,nx,1);
                uncertainty = reshape(uncertainty,nx,1);

                obj.nx= nx;
                obj.nu= nu;
                obj.ny= ny;
                obj.nw= nw;
                obj.M = [eye(obj.nx);-eye(obj.nx)];
                obj.m = [x0+uncertainty;uncertainty-x0];
                obj.k = 0;
                obj.nAd = 0;
                obj.nBd = 0;
                obj.nLd = 0;
                obj.nCd = 0;
                obj.nDd = 0;
                obj.nNd = 0;
                obj.Ad = zeros(obj.nx,0);
                obj.Bd = zeros(obj.nx,0);
                obj.Ld = zeros(obj.nx,0);
                obj.Cd = [];
                obj.Dd = [];
                obj.Nd = [];
                obj.flag_Delta = false;
                obj.vol = prod(2*uncertainty);
            end
        end
        
        function empty = update(obj,A,B,L,C,D,N,u,y,wp,wm)
            if any(size(A) ~= obj.nx)
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent matrix A dimensions.');
                throw(ME);
            elseif size(B,1) ~= obj.nx || size(B,2) ~= obj.nu
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent matrix B dimensions.');
                throw(ME);
            elseif size(C,1) ~= obj.ny || size(C,2) ~= obj.nx
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent matrix C dimensions.');
                throw(ME);
            elseif size(D,1) ~= obj.ny || size(D,2) ~= obj.nu
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent matrix D dimensions.');
                throw(ME);
            elseif size(L,1) ~= obj.nx || size(L,2) ~= obj.nw
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent matrix L dimensions.');
                throw(ME);
            elseif size(N,1) ~= obj.ny || size(N,2) ~= obj.nw
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent matrix N dimensions.');
                throw(ME);
            end
            try 
                u = reshape(u,[obj.nu,1]);
            catch ME0
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent input vector u dimensions.');
                ME = addCause(ME,ME0);
                throw(ME);
            end
            try
                y = reshape(y,[obj.ny,1]);
            catch ME0
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent output vector y dimensions.');
                ME = addCause(ME,ME0);
                throw(ME);
            end
            try
                wp= reshape(wp,[obj.nw,1]);
            catch ME0
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent vector dimensions.');
                ME = addCause(ME,ME0);
                throw(ME);
            end
            try
                wm= reshape(wm,[obj.nw,1]);
            catch ME0
                ME = MException('SVO:update:inconsistentDimensions',...
                    'Inconsistent vector dimensions.');
                ME = addCause(ME,ME0);
                throw(ME);
            end
            
            if ~obj.flag_Delta
                %prevents further attempts to add disturbances
                obj.flag_Delta = true;
            end
            tol = obj.numTol;
            wp = wp+tol;
            wm = wm-tol;
            if obj.nCd > 0 && obj.k == 0
                [Mx,mx,empty] = obj.projection(1:obj.nx);
                if empty
                    ME = MException('SVO:update:invalidInitialization',...
                        'The set of initial conditions is empty.');
                    throw(ME);
                end
                [MX,mX] = obj.cvxHull([Mx mx],[-Mx mx]);
            else
                MX = [];
                mX = [];
            end
            
            if any([obj.nCd obj.nLd obj.nNd] > 0)
                Mw = [eye(obj.nw);-eye(obj.nw)];
                mw = [wp;-wm];
                try 
                    [MW,mW] = obj.cvxHull([Mw mw],[-Mw mw]);
                catch ME1
                    ME2 = MException('SVO:update:invalidInput',...
                        'The given upper/lower bounds of the disturbances are not valid');
                    ME = addCause(ME1,ME2);
                    throw(ME);
                end
            else
                MW = [];
                mW = [];
            end    
            
            if any([obj.nCd obj.nDd obj.nBd] > 0)
                MU = [1;-1];
                mU = [1;1];
            else
                MU = [];
                mU = [];
            end
                
            
%                 keyboard
            if obj.k == 0
                  obj.M = [blkdiag(blkdiagn(MX,obj.nCd),obj.M), zeros(size(obj.M,1)+size(MX,1)*obj.nCd,obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                           obj.Cd,C,N,zeros(obj.ny,obj.nLd*obj.nw+obj.nBd),obj.Dd*blkdiagn(u,obj.nDd),obj.Nd,zeros(obj.ny,obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                           -obj.Cd,-C,-N,zeros(obj.ny,obj.nLd*obj.nw+obj.nBd),-obj.Dd*blkdiagn(u,obj.nDd),-obj.Nd,zeros(obj.ny,obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                           zeros(2*obj.nw,obj.nx*(1+obj.nCd)),[eye(obj.nw);-eye(obj.nw)],zeros(obj.nw*2,obj.nBd+obj.nLd*obj.nw+obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                           zeros(size(MU,1)*obj.nBd,obj.nx*(1+obj.nCd)+obj.nw),blkdiagn(MU,obj.nBd),zeros(size(MU,1)*obj.nBd,obj.nLd*obj.nw+obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                           zeros(size(MW,1)*obj.nLd,obj.nx*(1+obj.nCd)+obj.nw+obj.nBd),blkdiagn(MW,obj.nLd),zeros(size(MW,1)*obj.nLd,obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                           zeros(size(MU,1)*obj.nDd,obj.nx*(1+obj.nCd)+obj.nw+obj.nBd+obj.nw*obj.nLd),blkdiagn(MU,obj.nDd),zeros(size(MU,1)*obj.nDd,obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                           zeros(size(MW,1)*obj.nNd,obj.nx*(1+obj.nCd)+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd),blkdiagn(MW,obj.nNd),zeros(size(MW,1)*obj.nNd,obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));                          
                           zeros(size(MW,1)*obj.nCd,obj.nx*(1+obj.nCd)+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd+obj.nNd*obj.nw),blkdiagn(MW,obj.nCd),zeros(size(MW,1)*obj.nCd,obj.nCd*(1+obj.nBd+obj.nw*obj.nLd));
                           zeros(size(MU,1)*obj.nCd,obj.nx*(1+obj.nCd)+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd+obj.nNd*obj.nw+obj.nw*obj.nCd),blkdiagn(MU,obj.nCd),zeros(size(MU,1)*obj.nCd,obj.nCd*(obj.nBd+obj.nw*obj.nLd));
                           zeros(size(MW,1)*obj.nCd*obj.nLd,obj.nx*(1+obj.nCd)+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd+obj.nNd*obj.nw+obj.nw*obj.nCd+obj.nCd),blkdiagn(MW,obj.nLd*obj.nCd),zeros(size(MW,1)*obj.nCd*obj.nLd,obj.nCd*obj.nBd);
                           zeros(size(MU,1)*obj.nCd*obj.nBd,obj.nx*(1+obj.nCd)+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd+obj.nNd*obj.nw+(obj.nw+1+obj.nLd*obj.nw)*obj.nCd),blkdiagn(MU,obj.nCd*obj.nBd)];
                  obj.m = [repmat(mX,[obj.nCd,1]);
                           obj.m;
                           y-D*u+tol;
                           -y+D*u+tol;
                           wp+tol;
                           -wm+tol;
                           repmat(mU,[obj.nBd,1]);
                           repmat(mW,[obj.nLd,1]);
                           repmat(mU,[obj.nDd,1]);
                           repmat(mW,[obj.nNd,1]);
                           repmat(mW,[obj.nCd,1]);
                           repmat(mU,[obj.nCd,1]);
                           repmat(mW,[obj.nCd*obj.nLd,1]);
                           repmat(mU,[obj.nCd*obj.nBd,1])];
            else
                Mnew = [eye(obj.nx), zeros(obj.nx,obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                    -eye(obj.nx), zeros(obj.nx,obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                    C,N,zeros(obj.ny,obj.nBd+obj.nLd*obj.nw),obj.Dd*blkdiagn(u,obj.nDd),obj.Nd,zeros(obj.ny,obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                    -C,-N,zeros(obj.ny,obj.nBd+obj.nLd*obj.nw),-obj.Dd*blkdiagn(u,obj.nDd),-obj.Nd,zeros(obj.ny,obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                    zeros(2*obj.nw,obj.nx),[eye(obj.nw);-eye(obj.nw)],zeros(2*obj.nw,obj.nBd+obj.nLd*obj.nw+obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                    zeros(size(MU,1)*obj.nBd,obj.nx+obj.nw),blkdiagn(MU,obj.nBd),zeros(size(MU,1)*obj.nBd,obj.nLd*obj.nw+obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                    zeros(size(MW,1)*obj.nLd,obj.nx+obj.nw+obj.nBd),blkdiagn(MW,obj.nLd),zeros(size(MW,1)*obj.nLd,obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                    zeros(size(MU,1)*obj.nDd,obj.nx+obj.nw+obj.nBd+obj.nw*obj.nLd),blkdiagn(MU,obj.nDd),zeros(size(MU,1)*obj.nDd,obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                    zeros(size(MW,1)*obj.nNd,obj.nx+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd),blkdiagn(MW,obj.nNd),zeros(size(MW,1)*obj.nNd,obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd));
                    zeros(size(MW,1)*obj.nCd,obj.nx+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd+obj.nNd*obj.nw),blkdiagn(MW,obj.nCd),zeros(size(MW,1)*obj.nCd,obj.nCd*(1+obj.nBd+obj.nw*obj.nLd));
                    zeros(size(MU,1)*obj.nCd,obj.nx+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd+obj.nNd*obj.nw+obj.nw*obj.nCd),blkdiagn(MU,obj.nCd),zeros(size(MU,1)*obj.nCd,obj.nCd*(obj.nBd+obj.nw*obj.nLd));
                    zeros(size(MW,1)*obj.nCd*obj.nLd,obj.nx+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd+obj.nNd*obj.nw+obj.nw*obj.nCd+obj.nCd),blkdiagn(MW,obj.nLd*obj.nCd),zeros(size(MW,1)*obj.nCd*obj.nLd,obj.nCd*obj.nBd);
                    zeros(size(MU,1)*obj.nCd*obj.nBd,obj.nx+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd+obj.nNd*obj.nw+(obj.nw+1+obj.nLd*obj.nw)*obj.nCd),blkdiagn(MU,obj.nCd*obj.nBd);                  
                    zeros(size(obj.M,1),obj.nx+obj.nw+obj.nBd+obj.nw*obj.nLd+obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nw*obj.nLd))];
                Mold = [-obj.Ad,zeros(obj.nx,obj.nx*obj.nCd*(1+obj.nAd)),-obj.A,-obj.L,-obj.Bd*blkdiagn(u,obj.nBd),-obj.Ld,zeros(obj.nx,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw))*obj.k+(obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*(obj.k-1));
                    obj.Ad,zeros(obj.nx,obj.nx*obj.nCd*(1+obj.nAd)),obj.A,obj.L,obj.Bd*blkdiagn(u,obj.nBd),obj.Ld,zeros(obj.nx,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw))*obj.k+(obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*(obj.k-1));
                    zeros(obj.ny,obj.nx*obj.nAd),obj.Cd*blkdiagn(obj.Ad,obj.nCd),obj.Cd*blkdiagn(obj.A,obj.nCd),zeros(obj.ny,obj.nx+obj.nw+obj.nBd+obj.nDd+obj.nw*(obj.nLd+obj.nNd)),obj.Cd*blkdiagn(obj.L,obj.nCd),obj.Cd*blkdiagn(obj.B*obj.u,obj.nCd),obj.Cd*blkdiagn(obj.Ld,obj.nCd),obj.Cd*blkdiagn(obj.Bd*blkdiagn(obj.u,obj.nBd),obj.nCd),zeros(obj.ny,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*(obj.k-1));
                    zeros(obj.ny,obj.nx*obj.nAd),-obj.Cd*blkdiagn(obj.Ad,obj.nCd),-obj.Cd*blkdiagn(obj.A,obj.nCd),zeros(obj.ny,obj.nx+obj.nw+obj.nBd+obj.nDd+obj.nw*(obj.nLd+obj.nNd)),-obj.Cd*blkdiagn(obj.L,obj.nCd),-obj.Cd*blkdiagn(obj.B*obj.u,obj.nCd),-obj.Cd*blkdiagn(obj.Ld,obj.nCd),-obj.Cd*blkdiagn(obj.Bd*blkdiagn(obj.u,obj.nBd),obj.nCd),zeros(obj.ny,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*(obj.k-1));
                    zeros(2*obj.nw,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*obj.k);
                    zeros(size(MU,1)*obj.nBd,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*obj.k); 
                    zeros(size(MW,1)*obj.nLd,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*obj.k);
                    zeros(size(MU,1)*obj.nDd,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*obj.k);
                    zeros(size(MW,1)*obj.nNd,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*obj.k); 
                    zeros(size(MW,1)*obj.nCd,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*obj.k);
                    zeros(size(MU,1)*obj.nCd,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*obj.k);
                    zeros(size(MW,1)*obj.nCd*obj.nLd,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*obj.k);
                    zeros(size(MU,1)*obj.nCd*obj.nBd,(obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nAd*obj.nx+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw+obj.nCd*obj.nx*(1+obj.nAd))*obj.k);
                    obj.M];
                obj.M = [Mnew Mold];
                obj.m = [obj.B*obj.u+tol;
                    -obj.B*obj.u+tol;
                    y-D*u+tol;
                    -y+D*u+tol;
                    wp;
                    -wm;
                    repmat(mU,[obj.nBd,1]);
                    repmat(mW,[obj.nLd,1]);
                    repmat(mU,[obj.nDd,1]);
                    repmat(mW,[obj.nNd,1]);
                    repmat(mW,[obj.nCd,1]);
                    repmat(mU,[obj.nCd,1]);
                    repmat(mW,[obj.nCd*obj.nLd,1]);
                    repmat(mU,[obj.nCd*obj.nBd,1]);
                    obj.m];
            end
            if obj.k+1 > obj.n && strcmpi(obj.accelerator,'off')
                nidx = (obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw)*obj.k+...
                    (obj.nAd+obj.nCd*(1+obj.nAd))*obj.nx*(obj.k-1);
                [H,h,empty] = obj.projection(1:nidx);
                if ~empty
                    obj.M = H;
                    obj.m = h;
                else
                    obj.M = [];
                    obj.m = [];
                end
            elseif obj.k+1 > obj.n && strcmpi(obj.accelerator,'on')
                nidx = obj.nDd+obj.nw*obj.nNd+obj.nCd*(obj.nw+1+obj.nBd+obj.nLd*obj.nw)+obj.nx+obj.nw+obj.nBd+obj.nLd*obj.nw;
                [H,h,empty] = obj.projection(1:nidx);
                if ~empty
                    obj.M = H;
                    obj.m = h;
                else
                    obj.M = [];
                    obj.m = [];
                end
            else
                empty = obj.isempty();
            end
            if ~empty && (obj.nAd > 0 || (obj.nCd >0 && obj.k > 0) )
                %CAREFUL: this might generate an error because if the set
                %is empty, the number of columns in MX will not match the
                %number of uncertaintie
                [Mx,mx,void] = obj.projection(1:obj.nx);
                if void
                    warning('SVO:update:The SV estimate will be considered empty. It is possible that the solver ran into numerical problems. Consider increasing the tolerances.')
                    empty = 1;
                else
                    [MX,mX] = obj.cvxHull([Mx mx],[-Mx mx]);
                end
            else
                MX = [];
                mX = [];
            end
            obj.M = blkdiag(blkdiagn(MX,obj.nAd*(1+obj.nCd)+obj.nCd*(obj.k>0)),obj.M);
            obj.m = [repmat(mX,[obj.nAd*(1+obj.nCd)+obj.nCd*(obj.k>0) 1]);obj.m];
            
            if obj.k+1<= obj.n
                obj.k = obj.k + 1;
            elseif strcmpi(obj.accelerator,'on')
                obj.k = 1;
            end
            
            obj.A = A;
            obj.B = B;
            obj.L = L;
            obj.u = u;
        end
        
        function empty = isempty(obj)
            [~,~,status] = obj.LPsolve(zeros(1,size(obj.M,2)),obj.M,obj.m);
            if status == 1
                empty = false;
            else
                empty = true;
            end
        end
        
        function [Mx,mx,empty] = projection(obj,idx)
            if strcmpi(obj.approximation,'none')
                h = fourier([obj.M(:,1:obj.nx),obj.M(:,obj.nx+1:end),obj.m],...
                    idx,obj.numTol,obj.angTol);
                empty = obj.isempty();
            else
                [h,empty] = obj.boundingBox(obj.M,obj.m,idx);
            end
            if ~empty
                Mx = h(:,1:end-1);
                mx = h(:,end);  
            else
                Mx = [];
                mx = [];
            end
        end
        
        function addDelta(obj,varargin)
            if obj.flag_Delta
                ME = MException('SVO:addDelta:invalidAction',...
                    'Cannot add disturbances after a call to update has been made');
                throw(ME);
            elseif rem(numel(varargin),2) ~= 0
                ME = MException('SVO:addDelta:invalidInput',...
                    'The input must be a collection of pairs (id,mtx) with id=A,C and size(mtx)=size(id)');
                throw(ME)
            end
            for I = 1:2:numel(varargin)
                id = varargin{I};
                mtx= varargin{I+1};
                switch id
                    case 'A'
                        if any(size(mtx) ~= obj.nx)
                            ME = MException('SVO:addDelta:inconsistentDimensions',...
                                'The matrix A must be a %d x %d matrix',obj.nx,obj.nx);
                            throw(ME);
                        end
                        obj.nAd = obj.nAd + 1;
                        obj.Ad  = [obj.Ad,mtx];
                    case 'B'
                        if any(size(mtx) ~= [obj.nx obj.nu])
                            ME = MException('SVO:addDelta:inconsistentDimensions',...
                                'The matrix B must be a %d x %d matrix',obj.nx,obj.nu);
                            throw(ME);
                        end
                        obj.nBd = obj.nBd + 1;
                        obj.Bd  = [obj.Bd,mtx];
                    case 'C'
                        if any(size(mtx) ~= [obj.ny obj.nx])
                            ME = MException('SVO:addDelta:inconsistentDimensions',...
                                'The matrix C must be a %d x %d matrix',obj.ny,obj.nx);
                            throw(ME);
                        end
                        obj.nCd = obj.nCd + 1;
                        obj.Cd  = [obj.Cd,mtx];
                    case 'D'
                        if any(size(mtx) ~= [obj.ny obj.nu])
                            ME = MException('SVO:addDelta:inconsistentDimensions',...
                                'The matrix D must be a %d x %d matrix',obj.ny,obj.nu);
                            throw(ME);
                        end
                        obj.nDd = obj.nDd + 1;
                        obj.Dd  = [obj.Dd,mtx];
                    case 'L'
                        if any(size(mtx) ~= [obj.nx obj.nw])
                            ME = MException('SVO:addDelta:inconsistentDimensions',...
                                'The matrix L must be a %d x %d matrix',obj.nx,obj.nw);
                            throw(ME);
                        end
                        obj.nLd = obj.nLd + 1;
                        obj.Ld  = [obj.Ld,mtx];
                    case 'N'
                        if any(size(mtx) ~= [obj.nx obj.nw])
                            ME = MException('SVO:addDelta:inconsistentDimensions',...
                                'The matrix N must be a %d x %d matrix',obj.ny,obj.nw);
                            throw(ME);
                        end
                        obj.nNd = obj.nNd + 1;
                        obj.Nd  = [obj.Nd,mtx];
                    otherwise
                        ME = MException('SVO:addDelta:invalidInput',...
                            'The input must be a collection of pairs (id,mtx) with id=A,C,D and size(mtx)=size(id)');
                        throw(ME);
                end                        
            end
        end
        
        function reset(obj,x,u)
            if numel(x) ~= obj.nx
                ME = MException('SVO:reset:inconsistentDimensions',...
                    'The number of state variables nx does not match the number of elements in x.');
                throw(ME);
            elseif numel(u) ~= obj.nx
                ME = MException('SVO:reset:inconsistentDimensions',...
                    'The number of state variables nx does not match the number of elements in u.');
                throw(ME);
            end
            x = reshape(x,obj.nx,1);
            u = reshape(u,obj.nx,1);
            obj.k = 0;
            obj.M = [eye(obj.nx);-eye(obj.nx)];
            obj.m = [x+u;u-x];
        end
        function [x,u] = getStateEstimate(obj)
            if ~strcmpi(obj.approximation,'hypercube')
                ME = MException('SVO:getStateEstimate:invalidSVOsetting',...
                    'getStateEstimate requires the SVO to have the property "approximation" set to "hypercube"');
                throw(ME);
            end
            a = 1+obj.nAd*obj.nx*(obj.nCd+1)+obj.nCd*obj.nx;
            b = a+obj.nx-1;
            [~,h] = obj.projection(a:b);
            x_upr = h(1:obj.nx);
            x_lwr = -h(obj.nx+1:2*obj.nx);
            x = mean([x_upr,x_lwr],2);
            u = -diff([x_upr,x_lwr],1,2);  
        end
    end
end
            