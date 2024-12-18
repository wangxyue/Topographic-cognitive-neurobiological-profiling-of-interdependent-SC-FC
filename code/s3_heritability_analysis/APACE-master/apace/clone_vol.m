function V = clone_vol(Vtemplate, fname, descrip)
% Clone volume based on a template; create cloned volume
% FORMAT V = clone_vol(Vtemplate, fname, descrip)
%
% Clones a volume, while discarding any information we don't want to
% clone.
%_______________________________________________________________________
% Version: http://github.com/nicholst/APACE/tree/34ea883
%          2018-09-21 09:59:07 +0100

if iscell(Vtemplate)
    V = Vtemplate{1};
else
    V = Vtemplate(1);
end
V        = rmfield(V,{'fname','descrip','n','private'});
V.fname  = fname;
V.dim    = V.dim(1:3);
V.dt     = [spm_type('float32'), V.dt(2)];
V.pinfo  = [1 0 0]';
V.descrip= descrip;

V=spm_create_vol(V);

return

