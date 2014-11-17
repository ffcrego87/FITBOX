function [M,m] = reduce_rows(A,b)

zerotol = 1e-6;

%normalize
a = sqrt(sum(A.*A,2));
% A = A./a(:,ones(1,n));
% b = b./a;
%remove zero rows
idx = a>zerotol;
A = A(idx,:);
b = b(idx,:);
[m,n] = size(A);
aux = A*A';
%check if full-dimensional
[~,fval,exitflag,~,lambda] = linprog(-[zeros(1,n) 1],[A ones(m,1)],b);
if exitflag < 0
    ME = MException('reduce_rows:LP_error',...
        'LINGPROG failed to execute. Exitflag = %i',exitflag);
    throw(ME);
end

R=-fval;

if R < zerotol 
    %polytope is not full-dimensional
    E = lambda.ineqlin > zerotol;
    Ae= A(E,:);
    be= b(E);
else
    E = false(m,1);
    Ae = [];
    be = [];
end

idx = false(m,1) & ~E; 
for i = 1:m
    if ~E(i)
        jdx = abs(A(~E,:)*A(i,:)')>zerotol;
        [~,~,exitflag,~,lambda] = linprog(-1,A(jdx,:)*A(i,:)',b(jdx));
        if exitflag < 0
            ME = MException('reduce_rows:LP_error',...
                'LINGPROG failed to execute. Exitflag = %i',exitflag);
            throw(ME);
        end
        temp = false(m,1);
        temp(jdx) = lambda.ineqlin > zerotol;    
        idx = idx | temp;
    end
end

idx = idx | E;

M = A(idx,:);
m = b(idx,:);
