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
m1 = 0.001; m2 = 0.0004; m3 = 0.0019; m456=0.0004+0.0055+0.0065; %m4 = 0.0004; m5 = 0.0055; m6 = 0.0065;
I1 = 1; I2 = 1; I3 = 1; I4 = 1; I5 = 1;
k1=200; k2 = 375; k3 = 425.746*1000; k4 = 142.8571; %front, back,rotational, interdigital
cv11 = 1; cv12 = 1; cv13 = 1; cv2 = 1; cv3 = 1; cv4 = 1;

g = [ 0 , 0 , -9.81]; % gravity vector

%inputs:
%     cv11 cv12 cv13 cv2 cv31 cv32
F1=-9; 
d1 = 2; d2 = 2; d3 = 2; d4 = 4;
dc1=0.0225; dc2 = 0.0048; dc3 = 0.014055; dc4 = 0.010375; dc5 = 0.010325; dc6=0.011;
ds1 = 0.045; ds2 = 0.0096; ds3 = 0.02811; ds4 = 0.02075; ds5 = 0.02065; ds6 = 0.011; ds7 = 0.02289; ds8 = 0.0423;
x04 = 0.2; 

lc = [ 0 0 0 1 0 ;  0 0 -dc1 0 0 ; 0 dc2 0  0 0 ; 0 0 -dc3 0 0 ; 0 0 -ds4 1 0; 0 0 -ds5 1 0;0 0 -dc6 0 0; -x04 0 0 1 0 ; -ds8 0 0 2 0; ds7 0 0 2 7  ];%ext force(XY table),mass,mass,mass,ext force(antag. spring), 
                                                                           %ext force(ligament),mass456, friction force, constraint force at the front, constraint force at the back
m = [ m1 , m2 , m3 , m456];
I = sym ( zeros ( 3 , 3 , 4 ) );
I(:,:,1) = I1 * eye ( 3 ); I(:,:,2) = I2 * eye ( 3 ); I(:,:,3) = I3 * eye ( 3 ); I(:,:,4) = I4 * eye ( 3 ); %I(:,:,5) = I5 * eye ( 3 );
% j = sym ( zeros ( 3 , 5 , 3 ) );
j = sym ( zeros ( 1 , 5 , 10 ) );
% j(:,:,1) = [ 3 inf inf inf inf ; 1 inf 0 0 0 ; 2 inf 0 0 0 ];
j(1,:,1) = [ 0 -inf 0 0 0 ];%displacement due to the force
j(1,:,2) = [ 0 0 0 0 -ds1 ];%%%%%I do not know how to add the spring here
j(1,:,3) = [ 0 0 0 ds2 0 ];
j(1,:,4) = [ 1 inf 0 0 -ds3 ];%joint 3 movement around x axis
j(1,:,5) = [ 3 inf 0 0 -ds4 ];%joint 4 rotation around z
j(1,:,6) = [ 0 0 0 0 -ds5 ];
j(1,:,7) = [ 0 0 0 0 -ds6 ];
j(1,:,8) = [ 2 inf -x04 0 0];%rotation around y axis
j(1,:,9) = [ 0 0 -ds8+x04 0 0 ];
j(1,:,10) = [ 0 0 ds8+ds7 0 0 ];

% jkd = sym (zeros ( 3 , 2 , 9 ) );
jkd = sym (zeros ( 3 , 2 , 4 ) );
% jkd(2,:,4) = [ cv11 0 ]; jkd(2,:,5) = [ cv12 0 ]; jkd(2,:,6) = [ cv13 0 ];
% jkd(1:2,:,7) = [ k2 0;  cv2 0 ]; jkd(1:2,:,8) = [ k31 0; cv31 0 ]; jkd(1:2,:,9) = [ k32 0; cv32 0 ];
jkd(3,:,1) = [ F1 -9 ]; jkd(2:3,:,2) = [ cv2 0; 10 10 ]; jkd(1:2,:,3) = [ k3 0; cv3 0 ]; jkd(1:2,:,4) = [ k4 0; cv4 0 ];%the last one is a inetraction between the antagonistic springs

% EOM:h
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM5_S ( M , T , Dd , fg , fj , rcn , Tcn , Dcn , qf , uf , 1 );
plot ( t , z );
pause;

% animation
AnimEOM ( t , z , rj , qf , uf );