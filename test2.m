% EOM Simulation:

% properties:
lc1 = 1; lc2 = 1; lc3 = 1; lc4 = 1; % link COM pos.
lf1 = 1; lcn1 = 1;
m1 = 1; m2 = 1; m3 = 1; m4 = 1; % mass
I1 = 1e-3; I2 = 1e-3; I3 = 1e-3; I4 = 1e-3; % rotational inertia
l1 = 2; l2 = 2; % link length
k1 = 1e-1; k2 = 1e-1; k3 = 1e-1; k4 = 1e-1; % joint spring coeff.
g = [ 0 , 0 , -9.81]; % gravity vector

%inputs:
lc = [ lc1 0 0 0 0 ; lc2 0 0 0 0 ; lc3 0 0 0 0 ; lc4 0 0 0 1 ; lf1 0 0 1 3 ; lcn1 0 0 2 3 ];
m = [ m1 , m2 , m3 , m4 ];
I = sym ( zeros ( 3 , 3 , 4 ) );
I(:,:,1) = I1 * eye ( 3 ); I(:,:,2) = I2 * eye ( 3 ); I(:,:,3) = I3 * eye ( 3 ); I(:,:,4) = I4 * eye ( 3 );
j = sym ( zeros ( 1 , 5 , 3 ) );
j(:,:,1) = [ 2 inf 0 0 0 ]; j(:,:,2) = [ 2 inf l1 0 0 ]; j(:,:,3) = [ 2 inf l2 0 0 ];
j(:,:,4) = [ 2 inf l1 0 0 ]; % 2 joints coincides at 2nd joint (at the end of 1st link) with different rotation angle
j(:,:,5) = zeros ( 1 , 5 ); j(:,:,6) = zeros ( 1 , 5 ); % force and constraint exerted in 3rd link frame
jkd = sym (zeros ( 3 , 2 , 4 ) );
jkd(1,:,1) = [ k1 0 ]; jkd(1,:,2) = [ k2 0 ]; jkd(1,:,3) = [ k3 0 ]; jkd(1,:,3) = [ k4 0 ];

% EOM:
[ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g );

% numerical simulation
[ t , z , tfinal ] = SimEoM ( M , T , Dd , fg , fj , qf , uf , 7 );
plot ( t , z );
pause;

% animation
AnimEOM ( t , z , rj , qf , uf );

