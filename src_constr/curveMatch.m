function [crv_1,crv_2,u_order,u_knotvctr]=curveMatch(crv_1,crv_2)
% match curve degree and knotvctr
%
u_order_tar=max(crv_1.u_order,crv_2.u_order);
crv_1=crv_1.addOrder(u_order_tar-crv_1.u_order);
crv_2=crv_2.addOrder(u_order_tar-crv_2.u_order);

% normalize knotvctr
crv_1.u_knotvctr=(crv_1.u_knotvctr-crv_1.u_knotvctr(1))/(crv_1.u_knotvctr(end)-crv_1.u_knotvctr(1));
crv_2.u_knotvctr=(crv_2.u_knotvctr-crv_2.u_knotvctr(1))/(crv_2.u_knotvctr(end)-crv_2.u_knotvctr(1));

u_knotvctr_1=crv_1.u_knotvctr;
u_knotvctr_2=crv_2.u_knotvctr;

% merge the knot vectors of u
knots_com=unique([u_knotvctr_1,u_knotvctr_2]);
knots_1_ins=[];
knots_2_ins=[];
for i=1:length(knots_com)
    i1=sum(u_knotvctr_1 == knots_com(i));
    i2=sum(u_knotvctr_2 == knots_com(i));
    m=max(i1,i2);
    knots_1_ins=[knots_1_ins,knots_com(i)*ones(1,m-i1)];
    knots_2_ins=[knots_2_ins,knots_com(i)*ones(1,m-i2)];
end

if ~isempty(knots_1_ins),crv_1=crv_1.insertKnot(knots_1_ins);end
if ~isempty(knots_2_ins),crv_2=crv_2.insertKnot(knots_2_ins);end
u_order=crv_1.u_order;
u_knotvctr=crv_1.u_knotvctr;
end
