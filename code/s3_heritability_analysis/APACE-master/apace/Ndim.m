function nd = Ndim(X)
% Generic support function
% Returns the number of non-trivial dimensions
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/34ea883
%          2018-09-21 09:59:07 +0100

if isa(X,'nifti')
    nd = zeros(length(X),1);
    for i=1:length(X)
        nd(i) = X.hdr.dim(1);
    end
elseif isa(X,'gifti') && isfield(X,'cdata')
    nd = length(size(X.cdata));
    if ( nd==2 && size(X.cdata,2)==1 )
        nd = 1;
    end
elseif isnumeric(X)
    nd = length(size(X));
    if ( nd==2 && size(X,2)==1 )
        nd = 1;
    end
else
    error('Unknown image type')
end

return

