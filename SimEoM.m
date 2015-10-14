%% EOM Numerical Integrator:
% ====================================================

function [ t , z , tfinal ] = SimEoM ( M , T , D , fg , fj , q , u , dt )

global Mt;
global Tt;
global Dt;
global fgt;
global fjt;
global qt;
global ut;

Mt = M; Tt = T; Dt = D; fgt = fg; fjt = fj; qt = q; ut = u;

[ n , m ] = size ( T );
par = m;
z0 = 1e-6 * ones ( 1 , 2 * m );
t0 = 0;

% Standard ODE solver:
options = odeset (); %,'abstol',1*1e-6,'reltol',1*1e-6);
tspan = linspace ( t0 , t0 + dt , 500);
[ t , z , tfinal ] = ode113 ( @EOM , tspan , z0 , options , par );


function dz = EOM ( t , z , par )

global Mt;
global Tt;
global Dt;
global fgt;
global fjt;
global qt;
global ut;

qu = [ qt , ut ];
M = subs (  Mt , qu , z ); 
T = subs ( Tt , qu , z ); 
D = subs ( Dt , qu , z ); 
fg = subs ( fgt , qu , z ); 
fj = subs ( fjt , qu , z );

u = z( par + 1 : end );
dzt = double ( ( T.' * M * T ) \ ...
    ( T.' * ( fg - M * D * u ) + fj ) );
dz = [ u ; dzt ];

