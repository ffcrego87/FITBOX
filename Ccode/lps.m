function [v,x,status] = lps(c,a,b,e,objsen,lb,ub)
% LPS Linear program solver.
%
% [v,x,status] = lps(c,a,b,e,objsen,lb,ub)
%
% Interface to mex-file lpsmex
%
% Solves linear program
%
% {min,max} c^Tx  S/T  a*x ~ b
%
%  v := optimal cost     (NaN if infeasible, (+/-) inf if unbounded)
%  x := optimal solution (NaN if infeasible, (+/-) inf if unbounded)
%  status := solution status
%            1 if all OK, 2/3 if infeasible/unbounded
%  c,a,b := problem data (a can be sparse)
%  e := inequality determinator on constraints
%         less than: e(i) < 0 (default)
%         greater than: e(i) > 0
%         equals: e(i) = 0
%  objsen := min/max determinator
%               max: objsen < 0 (default)
%               min: objsen > 0
%  lb := lowerbounds on x (default = -inf)
%  ub := upper bounds on x (default = inf)
%
% Copyright 2015 Jeff Shamma
%  Address:
%  School of Electrical and Computer Engineering
%  Georgia Institute of Technology
%  777 Atlantic Dr NW
%  Atlanta, GA 30332-0250

[ncon,nvar] = size(a);

if nargin < 4,
   e = -ones(ncon,1);
end
if nargin < 5,
   objsen = -1;
end
if nargin < 6,
   lb = -ones(nvar,1)*inf;
end
if nargin < 7,
   ub = ones(nvar,1)*inf;
end

%
% Check dimensions
%

if length(c) ~= nvar,
   error('Dimension mismatch in a or c');
end
if length(b) ~= ncon,
   error('Dimension mismatch in a or b');
end
if length(e) ~= ncon,
   error('Dimension mismatch in a or e');
end
if length(lb) ~= nvar,
   error('Dimension mismatch in a or lb');
end
if length(ub) ~= nvar,
   error('Dimension mismatch in a or ub');
end

%
% Convert "a" to cplex sparse data structure
%

[matbeg,matcnt,matind,matval] = cpxprep(a);
numnz = length(matval);


%
% Call cplex interface
%

c = full(c(:)); % Mex-file does not use sparse structures
b = full(b(:));
e = full(e);
lb = full(lb);
ub = full(ub);

% Gurda Arslan added on 7/31/2002
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
matbeg = full(matbeg);
matcnt = full(matcnt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[v,x,status] = lpsmex(c,matbeg,matcnt,matind,matval,...
                          b,e,objsen,lb,ub,ncon,nvar,numnz);

if (v+1) == (v-1),  % CPX infinity
  v = inf*sign(v);
end

if status == 11,
  status = 1;
elseif status == 2,
  v = NaN;
  x = x*NaN;
elseif status == 3,
  v = inf*sign(-objsen);
  x = ones(nvar,1)*inf*sign(-objsen);
end





