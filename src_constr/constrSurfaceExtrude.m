function srf=constrSurfaceExtrude(crv,vctr)
% construct a surface by extruding a curve
%
coefs_v1=reshape(crv.coefs,crv.u_coef_num,1,crv.coef_dim);
crv=crv.translate(vctr);
coefs_v2=reshape(crv.coefs,crv.u_coef_num,1,crv.coef_dim);
coefs=cat(2,coefs_v1,coefs_v2);
srf=Surface(coefs,{crv.u_knotvctr,[0,0,1,1]});
end

