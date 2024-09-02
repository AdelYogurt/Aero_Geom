function vol=constrVolumeExtrud(srf,vctr)
% construct a volume by extruding a surface
%
coefs_v1=reshape(srf.coefs,srf.u_coef_num,srf.v_coef_num,1,srf.coef_dim);
srf=srf.translate(vctr);
coefs_v2=reshape(srf.coefs,srf.u_coef_num,srf.v_coef_num,1,srf.coef_dim);
coefs=cat(3,coefs_v1,coefs_v2);
vol=Volume(coefs,{srf.u_knotvctr,srf.v_knotvctr,[0,0,1,1]});
end

