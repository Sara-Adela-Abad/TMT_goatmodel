% EOM Simulation:
clc
clear all
close all
pause(1e-2)

% properties:
% syms d1 d2 d3...
%     dc2 dc3...
%     x03 xc3...
%     m1 m2 m3...
%     I1 I2 I3...
%     k2 k31 k32...
%     cv11 cv12 cv13 cv2 cv31 cv32
F1=-9; 
d1 = 2; d2 = 2; d3 = 2; d4 = 4;
dc1=0.5; dc2 = 0.5; dc3 = 0.5; dc4 = 3;
ds1 = 1; ds2 = 1; ds3 = 1;
x04 = 0.2; xc4 = 0.2;
m1 = 1; m2 = 1; m3 = 1; m4 = 1; m5 = 1;
I1 = 1; I2 = 1; I3 = 1; I4 = 1; I5 = 1;
k2 = 1; k3 = 1; k4 = 1;
cv11 = 1; cv12 = 1; cv13 = 1; cv2 = 1; cv3 = 1; cv4 = 1;

g = [ 0 , 0 , -9.81]; % gravity vector

%inputs:
lc = [ 0 0 0 0 0 ;  0 0 -dc1 0 0 ; 0 dc2 0  0 0 ; 0 -dc3  0 0 0 ; -dc4 0 0 0 0 ; x04-xc4 0 0 1 5 ; 0 0 -ds1 2 3 ; -ds2 0 0 2 5 ; ds3 0 0 2 5 ];
m = [ m1 , m2 , m3 , m4, m5];
I = sym ( zeros ( 3 , 3 , 5 ) );
I(:,:,1) = I1 * eye ( 3 ); I(:,:,2) = I2 * eye ( 3 ); I(:,:,3) = I3 * eye ( 3 ); I(:,:,4) = I4 * eye ( 3 ); I(:,:,5) = I5 * eye ( 3 );
% j = sym ( zeros ( 3 , 5 , 3 ) );
j = sym ( zeros ( 1 , 5 , 9 ) );
% j(:,:,1) = [ 3 inf inf inf inf ; 1 inf 0 0 0 ; 2 inf 0 0 0 ];
j(1,:,1) = [ 0 0 inf 0 0 ];
j(1,:,2) = [ 0 0 0 0 -d1 ];
j(1,:,3) = [ 1 inf 0 d2 0 ];
j(1,:,4) = [ 3 inf 0 0 -d3 ];
j(1,:,5) = [2 inf 0 0 -d4 ];
% jkd = sym (zeros ( 3 , 2 , 9 ) );
jkd = sym (zeros ( 3 , 2 , 4 ) );
% jkd(2,:,4) = [ cv11 0 ]; jkd(2,:,5) = [ cv12 0 ]; jkd(2,:,6) = [ cv13 0 ];
% jkd(1:2,:,7) = [ k2 0;  cv2 0 ]; jkd(1:2,:,8) = [ k31 0; cv31 0 ]; jkd(1:2,:,9) = [ k32 0; cv32 0 ];
jkd(3,:,1) = [ F1 -9 ]; jkd(2:3,:,3) = [ cv3 0; 10 10 ]; jkd(1:2,:,4) = [ k4 0; cv4 0 ];

% EOM:h
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM5_S ( M , T , Dd , fg , fj , rcn , Tcn , Dcn , qf , uf , 1 );
plot ( t , z );
pause;

% animation
AnimEOM ( t , z , rj , qf , uf );