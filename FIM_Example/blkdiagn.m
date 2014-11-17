function out = blkdiagn(X,n)
    if n == 0
        out = [];
    elseif n > 0
        out = blkdiag(X,blkdiagn(X,n-1));
    else
        error('Error::Blkdiagn')
    end
end