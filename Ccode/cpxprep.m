function [matbeg,matcnt,matind,matval] = cpxprep(a)
% CPXPREP Prepare matrix data for CPLEX.
%
% [matbeg,matcnt,matind,matval] = cpxprep(a)
%
% Copyright 2015 Jeff Shamma
%  Address:
%  School of Electrical and Computer Engineering
%  Georgia Institute of Technology
%  777 Atlantic Dr NW
%  Atlanta, GA 30332-0250

[ncon,nvar] = size(a);
a = sparse(a);

if nzmax(a) <= 1,     % Trick in case of a = 0
  holdit = a(1,1);
  a(1) = 1;
  fakeit = 1;
else
  fakeit = 0;
end

[matind,j,matval] = find(a);
numnz = length(matval);
matind = matind(:) - 1;
matind = [matind;zeros(nvar-numnz,1)];

% Edited by Paulo Rosa, 2009
matcnt = full(sum(spones(a), 1));
% End of editing

matbeg = cumsum(matcnt);
matbeg = [0 matbeg(1:nvar-1)];

if fakeit,
  matval(1) = holdit;
end
