function wc = wc(z)
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

%WC
%    WC = WC(Q3,Q4,Q5,U3,U4,U5)

%    This function was generated by the Symbolic Math Toolbox version 5.9.
%    14-Oct-2015 12:19:46

t2 = cos(q3);
t3 = sin(q3);
t4 = sqrt(2.0);
t5 = t4.*(1.0./4.0);
t6 = sqrt(6.0);
t8 = t6.*(1.0./4.0);
t7 = t5-t8;
t9 = t2.^2;
t10 = t5+t8;
t11 = t3.^2;
t12 = t7.*t9.*u3;
t13 = t7.*t11.*u3;
t14 = t12+t13;
t15 = t9.*t10.*u3;
t16 = t10.*t11.*u3;
t17 = t15+t16;
t18 = cos(q4);
t19 = sin(q4);
t20 = t10.*t19;
t27 = t3.*t7.*t18;
t21 = t20-t27;
t22 = t7.*t19;
t23 = t3.*t10.*t18;
t24 = t22+t23;
t25 = t7.*t18;
t31 = t3.*t10.*t19;
t26 = t25-t31;
t28 = t10.*t18;
t29 = t3.*t7.*t19;
t30 = t28+t29;
t32 = t30.*u4;
t57 = t2.*t7.*t18.*u3;
t33 = t32-t57;
t34 = t2.*t18.*t33;
t35 = t21.*u4;
t58 = t2.*t7.*t19.*u3;
t36 = t35-t58;
t37 = t2.*t19.*t36;
t38 = -t13+t34+t37;
t39 = t3.*t18.*u3;
t40 = t2.*t19.*u4;
t41 = t39+t40;
t42 = t24.*t41;
t43 = t2.*t18.*u4;
t59 = t3.*t19.*u3;
t44 = t43-t59;
t45 = t26.*t44;
t46 = t15+t42+t45;
t47 = t26.*u4;
t48 = t2.*t10.*t18.*u3;
t49 = t47+t48;
t50 = t21.*t49;
t51 = t24.*u4;
t52 = t2.*t10.*t19.*u3;
t53 = t51+t52;
t54 = t2.*t3.*t7.*t10.*u3;
t56 = t30.*t53;
t55 = t50+t54-t56;
t60 = sin(q5);
t61 = cos(q5);
t62 = t3.*t61;
t63 = t2.*t19.*t60;
t64 = t62+t63;
t65 = t3.*t60;
t78 = t2.*t19.*t61;
t66 = t65-t78;
t67 = t30.*t61;
t68 = t2.*t7.*t60;
t69 = t67+t68;
t70 = t26.*t60;
t71 = t2.*t10.*t61;
t72 = t70+t71;
t73 = t30.*t60;
t82 = t2.*t7.*t61;
t74 = t73-t82;
t75 = t26.*t61;
t77 = t2.*t10.*t60;
t76 = t75-t77;
t79 = t3.*t7.*t60;
t122 = t2.*t7.*t19.*t61;
t80 = t79-t122;
t81 = t80.*u3;
t83 = t74.*u5;
t84 = t21.*t61.*u4;
t85 = t81+t83+t84;
t86 = t66.*t85;
t87 = t3.*t7.*t61;
t88 = t2.*t7.*t19.*t60;
t89 = t87+t88;
t90 = t89.*u3;
t91 = t69.*u5;
t123 = t21.*t60.*u4;
t92 = t90+t91-t123;
t93 = t64.*t92;
t94 = -t34+t86+t93;
t95 = t2.*t60;
t96 = t3.*t19.*t61;
t97 = t95+t96;
t98 = t97.*u3;
t99 = t64.*u5;
t124 = t2.*t18.*t61.*u4;
t100 = t98+t99-t124;
t101 = t2.*t61;
t126 = t3.*t19.*t60;
t102 = t101-t126;
t103 = t102.*u3;
t104 = t2.*t18.*t60.*u4;
t127 = t66.*u5;
t105 = t103+t104-t127;
t106 = t72.*t105;
t125 = t76.*t100;
t107 = t42+t106-t125;
t108 = t3.*t10.*t60;
t128 = t2.*t10.*t19.*t61;
t109 = t108-t128;
t110 = t72.*u5;
t111 = t24.*t61.*u4;
t129 = t109.*u3;
t112 = t110+t111-t129;
t113 = t69.*t112;
t114 = t3.*t10.*t61;
t115 = t2.*t10.*t19.*t60;
t116 = t114+t115;
t117 = t116.*u3;
t118 = t24.*t60.*u4;
t130 = t76.*u5;
t119 = t117+t118-t130;
t120 = t74.*t119;
t121 = -t50+t113+t120;
wc = reshape([0.0,0.0,0.0,-t3.*t7.*t17+t3.*t10.*t14,-t24.*t38+t21.*t46-t2.*t18.*t55,t24.*t94+t21.*t107+t2.*t18.*t121,0.0,0.0,0.0,t7.*t14+t10.*t17,-t26.*t38+t30.*t46+t2.*t19.*t55,t76.*t94+t69.*t107+t66.*t121,0.0,0.0,0.0,t2.*t7.*t17-t2.*t10.*t14,-t3.*t55+t2.*t10.*t38+t2.*t7.*t46,-t72.*t94-t74.*t107+t64.*t121],[6, 3]);
