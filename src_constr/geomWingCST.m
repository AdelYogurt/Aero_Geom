function [srf_L,srf_U,srf_R,srf_T,crv_WR_L,crv_WR_U,crv_WT_L,crv_WT_U]=geomWingCST(CR,CT,SL,W,TU,TL,MU,ML)
% construct curve of boundary on wing
%
crv_WR_U=CurveCST(MU,[CR,TU]);
crv_WR_U=crv_WR_U.convertSpline();
crv_WR_U=crv_WR_U.insertDimension(1,0);

crv_WR_L=CurveCST(ML,[CR,-TL]);
crv_WR_L=crv_WR_L.convertSpline();
crv_WR_L=crv_WR_L.insertDimension(1,0);

crv_WT_U=CurveCST(MU,[CT,TU/CR*CT]);
crv_WT_U=crv_WT_U.convertSpline();
crv_WT_U=crv_WT_U.insertDimension(1,0);
crv_WT_U=crv_WT_U.translate([W*tan(SL),W,0]);

crv_WT_L=CurveCST(ML,[CT,-TL/CR*CT]);
crv_WT_L=crv_WT_L.convertSpline();
crv_WT_L=crv_WT_L.insertDimension(1,0);
crv_WT_L=crv_WT_L.translate([W*tan(SL),W,0]);

srf_L=constrSurfaceBound(crv_WR_L,crv_WT_L);
srf_U=constrSurfaceBound(crv_WR_U,crv_WT_U);
srf_R=constrSurfaceBound(crv_WR_L,crv_WR_U);
srf_T=constrSurfaceBound(crv_WT_L,crv_WT_U);
end
