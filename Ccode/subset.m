function  v = subset(in,S2)
%SUBSET Check if vector/polytope object is subset of another polytope object.
%
% v = subset(x,S)
%     -or-
% v = subset(S1,S2)
%
% Checks if vector x or polytope object S1 is contained in S2 by computing
%
%  v(i) =  M2(i,:)*x - m2(i)  
%
%      -or-
%
%  v(i) = max M2(i,:)*x - m2(i)  s/t  M1*x <= m1
%          x
%
% where [M1,m1] = get(S1) and [M2,m2] = get(S2).
%
%  v := vector of costs
%  x or S1, S2 := input vector/polytope objects
%
% Copyright 2015 Jeff Shamma
%  Address:
%  School of Electrical and Computer Engineering
%  Georgia Institute of Technology
%  777 Atlantic Dr NW
%  Atlanta, GA 30332-0250

[M2,m2] = get(S2);
M2 = full(M2);     % Needs full matrix

if isa(in,'polytope')

  [M1,m1] = get(in);

  if size(M1,2) ~= size(M2,2),
    error('Polytope dimension mismatch.')
  end

  [matbeg,matcnt,matind,matval] = cpxprep(M1);
  [ncon,nvar] = size(M1);
  numnz = length(matval);

  [v,stats] = subsetmex(M2,matbeg,matcnt,matind,matval,...
                  m1,ncon,nvar,numnz,length(m2));

  v = v(:) - m2;

else

  x = in;

  if length(x) ~= size(M2,2),
    error('Polytope/vector dimension mismatch.')
  end

  v = M2*x - m2;

end

