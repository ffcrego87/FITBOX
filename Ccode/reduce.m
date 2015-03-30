function Snew = reduce(S,tol)
% REDUCE Eliminate redundant constraints in polytope object.
%
% Snew = reduce(S,tol)
%
% Removes redundant contraints in inequality description
% M*x <= m for polytope object S. Sensitive to scaling of S.
%
%  Snew := output polytope object
%  S := input polytope object
%  tol := tolerance to remove redundancies (default 1.0e-6)
%
% Copyright 2015 Jeff Shamma
%  Address:
%  School of Electrical and Computer Engineering
%  Georgia Institute of Technology
%  777 Atlantic Dr NW
%  Atlanta, GA 30332-0250

if nargin < 2,
   tol = 1.0e-6;
end 

[M,m] = get(S);
M = full(M);             % Needs full matrix

if norm(M) < tol^2,      % Check for all zeros
  if min(m) > tol,
    Snew = polytope(zeros(1,size(M,2)),1);
    return
  end
end

% M = [M(1,:);M];  % CPLEX Bug? don't ask...
% m = [m(1);m];

[matbeg,matcnt,matind,matval] = cpxprep(M);
[ncon,nvar] = size(M);
numnz = length(matval);

[costs,stats,ikeep] = reducemex(M,matbeg,matcnt,matind,matval,...
                  m,ncon,nvar,numnz,tol);

if length(find(stats==2)) ~= 0,
  Snew = polytope;
  return
end

ikeep = find(ikeep > 0.5);
M = M(ikeep,:);
m = m(ikeep,:);

Snew = polytope(M,m);



