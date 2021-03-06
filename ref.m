function ref = ref(z)
q1 = z(1);
q2 = z(2);
q3 = z(3);
q4 = z(4);
q5 = z(5);
u1 = z(6);
u2 = z(7);
u3 = z(8);
u4 = z(9);
u5 = z(10);

%REF
%    REF = REF(Q1,Q2,Q3,Q4,Q5)

%    This function was generated by the Symbolic Math Toolbox version 5.9.
%    14-Oct-2015 12:19:41

t2 = sin(q3);
t3 = cos(q3);
t4 = sqrt(2.0);
t5 = t4.*(1.0./4.0);
t6 = sqrt(6.0);
t9 = t6.*(1.0./4.0);
t7 = t5-t9;
t8 = sin(q5);
t10 = sin(q4);
t11 = cos(q5);
t12 = q2-2.0./2.5e1;
t13 = t5+t9;
t14 = t6.*2.325e-2;
t15 = cos(q4);
t16 = t7.*t12;
t17 = t12.*t13;
t18 = t6.*(3.7e1./1.0e3);
t19 = t13.*t15;
t20 = t2.*t7.*t10;
t21 = t19+t20;
t22 = t7.*t15;
t27 = t2.*t10.*t13;
t23 = t22-t27;
t24 = t3.*t10.*(9.0./5.0e2);
t25 = t7.*t15.*(9.0./5.0e2);
t26 = t3.*t10.*t11.*(2.3e1./1.0e3);
t28 = t11.*t23.*(2.3e1./1.0e3);
ref = reshape([t2.*(-9.19e-2)-t2.*t11.*(1.3e1./1.0e3)-t3.*t8.*t10.*(1.3e1./1.0e3)+1.0./5.0e1,0.0,t2.*(-2.7e1./1.0e3)+1.0./5.0e1,t2.*(-5.51e-2)+t24+1.0./5.0e1,t2.*(-9.19e-2)+t26-t2.*t8.*(2.3e1./1.0e3)+1.0./5.0e1,t2.*(-5.51e-2)-t24+1.0./5.0e1,t2.*(-9.19e-2)-t26+t2.*t8.*(2.3e1./1.0e3)+1.0./5.0e1,q1-t4.*2.325e-2+t14+t16-t3.*t7.*9.19e-2+t8.*t21.*(1.3e1./1.0e3)-t3.*t7.*t11.*(1.3e1./1.0e3),q1-t4.*(3.7e1./1.0e3)+t16+t18,q1-t4.*2.325e-2+t14+t16-t3.*t7.*(2.7e1./1.0e3),q1-t4.*2.325e-2+t14+t16-t3.*t7.*5.51e-2-t13.*t15.*(9.0./5.0e2)-t2.*t7.*t10.*(9.0./5.0e2),q1-t4.*2.325e-2+t14+t16-t3.*t7.*9.19e-2-t11.*t21.*(2.3e1./1.0e3)-t3.*t7.*t8.*(2.3e1./1.0e3),q1-t4.*2.325e-2+t14+t16-t3.*t7.*5.51e-2+t13.*t15.*(9.0./5.0e2)+t2.*t7.*t10.*(9.0./5.0e2),q1-t4.*2.325e-2+t14+t16-t3.*t7.*9.19e-2+t11.*t21.*(2.3e1./1.0e3)+t3.*t7.*t8.*(2.3e1./1.0e3),t4.*(-2.325e-2)-t14+t17-t3.*t13.*9.19e-2-t8.*t23.*(1.3e1./1.0e3)-t3.*t11.*t13.*(1.3e1./1.0e3)-9.0./5.0e1,t4.*(-3.7e1./1.0e3)+t17-t18-9.0./5.0e1,t4.*(-2.325e-2)-t14+t17-t3.*t13.*(2.7e1./1.0e3)-9.0./5.0e1,t4.*(-2.325e-2)-t14+t17+t25-t3.*t13.*5.51e-2-t2.*t10.*t13.*(9.0./5.0e2)-9.0./5.0e1,t4.*(-2.325e-2)-t14+t17+t28-t3.*t13.*9.19e-2-t3.*t8.*t13.*(2.3e1./1.0e3)-9.0./5.0e1,t4.*(-2.325e-2)-t14+t17-t25-t3.*t13.*5.51e-2+t2.*t10.*t13.*(9.0./5.0e2)-9.0./5.0e1,t4.*(-2.325e-2)-t14+t17-t28-t3.*t13.*9.19e-2+t3.*t8.*t13.*(2.3e1./1.0e3)-9.0./5.0e1],[7, 3]);
