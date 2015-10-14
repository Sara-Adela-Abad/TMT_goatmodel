%% EOM Numerical Integrator:
% ====================================================

function [ t , z , tfinal ] = SimEoM5_S ( M , T , D , fg , fj , rcn , Tcn , Dcn , q , u , dt )

global Mt;
global Tt;
global Dt;
global Tcnt;
% global Dcnt;
global rcnt;
global fgt;
global fjt;
global qt;
global ut;

Mt = M; Tt = T; Dt = D; Tcnt = Tcn; fgt = fg; fjt = fj; qt = q; ut = u; rcnt = rcn;
% Dcnt = Dcn; 

[ n , m ] = size ( T );
par = m;
z0 = 0.1 * ones ( 1 , 2 * m );
% z0(1:m) = [0 0 0];
% z0(1) = 0;
% z0(m+1:m+3) = [0 0 0];
% z0(m+1) = 0;
t0 = 0;

% Standard ODE solver:
options = odeset (); %,'abstol',1*1e-6,'reltol',1*1e-6);
tspan = linspace ( t0 , t0 + dt , 500);
[ t , z , tfinal ] = ode113 ( @EOM , tspan , z0 , options , par );


function dz = EOM ( t , z , par )

global Mt;
global Tt;
global Dt;
global Tcnt;
% global Dcnt;
global rcnt;
global fgt;
global fjt;
global qt;
global ut;

qu = [ qt , ut ];
M = subs (  Mt , qu , z ); 
T = subs ( Tt , qu , z ); 
D = subs ( Dt , qu , z );
Tcn = subs ( Tcnt , qu , z ); 
% Dcn = subs ( Dcnt , qu , z ); 
fg = subs ( fgt , qu , z ); 
fj = subs ( fjt , qu , z );
rcn = subs ( rcnt , qu , z );

rs1 = rcn(1,:)-rcn(2,:);
rs2 = rcn(1,:)-rcn(3,:);
k1 = 1e0;
k2 = 1e0;

u = z( par + 1 : end );
dzt = double ( ( T.' * M * T ) \ ...
    ( T.' * ( fg - M * D * u ) + fj +...
    k1*(sqrt(rs1*rs1'))*((Tcn(:,:,1)-Tcn(:,:,2))'*[1; 1; 1])+...
    k2*(sqrt(rs2*rs2'))*((Tcn(:,:,1)-Tcn(:,:,3))'*[1; 1; 1])));
dz = [ u ; dzt ];

