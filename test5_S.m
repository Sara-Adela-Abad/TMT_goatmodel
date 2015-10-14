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

d1 = 1; d2 = 1; d3 = 1;
dc1= 0.2; dc2 = 0.5; dc3 = 0.5; dc4 = 0.5; dc5 = 0.5; dc6 = 0.5;%in meters?
x06 = 0.2; xc6 = 0.2;
m2 = 1; m3 = 1; m4 = 1; m6 = 1;
I1 = 1; I2 = 1; I3 = 1; I4 = 1; I6 = 1;
k2 = 1; k31 = 1; k32 = 1;
cv11 = 1; cv12 = 1; cv13 = 1; cv2 = 1; cv31 = 1; cv32 = 1;

g = [ 0 , 0 , -9.81]; % gravity vector

%inputs:
lc = [ 0 0 0 1 0 ; 0 -dc1 0 0 0 ; 0 0 -dc2 0 0 ; 0 -dc3 0 0 0 ; 0 -dc4 0 1 0; 0 -dc5 0 0 0 ; -x06 0 0 1 0 ; xc6 0 0 1 6 ; dc6 0 0 1 6 ];
m = [ m2 , m3, m4, m6 ];
I = sym ( zeros ( 3 , 3 , 4 ) );
I(:,:,1) = I1 * eye ( 3 ); I(:,:,2) = I2 * eye ( 3 ); I(:,:,3) = I3 * eye ( 3 ); I(:,:,4) = I4 * eye ( 3 );
% j = sym ( zeros ( 3 , 5 , 3 ) );
j = sym ( zeros ( 1 , 5 , 9 ) );
% j(:,:,1) = [ 3 inf inf inf inf ; 1 inf 0 0 0 ; 2 inf 0 0 0 ];
j(1,:,1) = [ 0 0 inf 0 0];
j(1,:,2) = [ 0 0 0 -dc1 0 ];
j(1,:,3) = [ 1 inf 0 0 -dc2 ]; 
j(1,:,4) = [ 2 inf 0 -dc3 0 ];
j(1,:,5) = [ 0 0 0 -dc4 0 ];
j(1,:,6) = [ 0 0 0 -dc5 0];
j(1,:,7) = [ 3 inf -x06 0 0];
j(1,:,8) = [ 3 inf xc6 0 0];
j(1,:,9) = [ 3 inf dc6 0 0];


%j(1:2,:,3) = [ 3 inf 0 0 -d2 ; 1 inf 0 0 0 ];
% jkd = sym (zeros ( 3 , 2 , 9 ) );
jkd = sym (zeros ( 3 , 2 , 6 ) );
% jkd(2,:,4) = [ cv11 0 ]; jkd(2,:,5) = [ cv12 0 ]; jkd(2,:,6) = [ cv13 0 ];
% jkd(1:2,:,7) = [ k2 0; cv2 0 ]; jkd(1:2,:,8) = [ k31 0; cv31 0 ]; jkd(1:2,:,9) = [ k32 0; cv32 0 ];
jkd(2,:,2) = [ cv2 0 ]; jkd(2,:,3) = [cv31 0 ]; jkd(1:2,:,4) = [ k32 0; cv32 0 ]; jkd(1:2,:,5) = [ k32 0; cv32 0 ]; jkd(1:2,:,6) = [ k32 0; cv32 0 ];

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM5_S ( M , T , Dd , fg , fj , qf , uf , 1 );
plot ( t , z );
pause;

% animation
AnimEOM ( t , z , rj , qf , uf );