%% EOM Numerical Integrator:
% ====================================================

function [ t , z , tfinal ] = SimEoM5_S_v0_3( M , T , Dd , fg , fj , ref , Tef , qf , uf , dt )

global Mt ;
global Tt ;
global Dt ;
global Teft ;
global reft ;
global fgt ;
global fjt ;
global qt ;
global ut ;

Mt = M ; Tt = T ; Dt = Dd ; Teft = Tef ; fgt = fg ; fjt = fj ; qt = qf ; ut = uf ; reft = ref ;

[ n , m ] = size ( T ) ;
par = m ;
z0 = 1e-3 * ones ( 1 , 2 * m ) ;
t0 = 0 ;

% Standard ODE solver:
options = odeset () ;  %,'abstol',1*1e-6,'reltol',1*1e-6) ;
tspan = linspace ( t0 , t0 + dt , 10) ;
[ t , z , tfinal ] = ode113( @EOM , tspan , z0 , options , par ) ;


function dz = EOM ( t , z , par )

t

global Mt ;
global Tt ;
global Dt ;
global Teft ;
global reft ;
global fgt ;
global fjt ;
global qt ;
global ut ;

k_s2 = 375; k_s3 = 200 ; k_s4 = 250 ;
k_g = 300 ; c_vg = 1 ;%for asphalt
% k_g = 1e2 ; c_vg = 1e2 ;
% k_g = 0 ; c_vg = 0 ;

%%friction force coefficients
meu = 0.85;%frition coefficient for rubber -dry concrete
meu_s= 1.7; %static meu
v_0= 4e-5 ;%stribeck speed
v_c= 1; %viscosity coefficient

l_s4 = 3e-2 ; l_s2 = 3e-2 ; l_s3 = 3.5e-2 ;
h = [ 0 0 -46.5e-2 ] ;
input = -1e-3 ;

u = z( par + 1 : end ) ;

qu = [ qt , ut ].' ;
M = double ( subs(  Mt , qu , z ) ) ; 
T = double ( subs( Tt , qu , z ) ) ; 
D = double ( subs( Dt , qu , z ) ) ;
fg = double ( subs( fgt , qu , z ) ) ;
fj = double ( subs( fjt , qu , z ) ) ;
Tef = double ( subs( Teft , qu , z ) ) ; 
ref = double ( subs( reft , qu , z ) ) ;

fj(1) = fj(1) + input ;


r_sg = [ 0 0 ref( 1 , 3 ) ] - h ;
r_sgv = sqrt( r_sg * r_sg' ) ;
v_g = Tef( : , : , 1 ) * u ;

r_s4 = ref( 2 , : ) - ref( 3 , : ) ;
r_s2 = ref( 4 , : ) - ref( 5 , : ) ;
r_s3 = ref( 6 , : ) - ref( 7 , : ) ;
r_sv4 = sqrt( r_s4 * r_s4' ) ;
r_sv2 = sqrt( r_s2 * r_s2' ) ;
r_sv3 = sqrt( r_s3 * r_s3' ) ;

NormalForce = k_g * r_sgv - c_vg * v_g(3);
Fra=NormalForce*sign( v_g(3))*(meu+(meu_s-meu)*exp(-(v_g(3)/v_0)^2));
Frb=v_c*NormalForce*v_g(3);
FrictionForce = Fra+Frb;
%  FrictionForce = - ( k_g * r_sgv - c_vg * v_g(3) ) * meu * sign( v_g(3) )
% ( k_g * r_sgv - c_vg * v_g(3) ) * ( [ 0 0 1 ] * Tef( : , : , 1 ) )'
%  k_g * r_sgv - c_vg * v_g(3) ) * meu * sign( v_g(3) ) * ( [ 0 -1 0 ] * Tef( : , : , 1 ) )'
% k_s4 * ( 1 - l_s4 / r_sv4 ) * ( r_s4 * ( Tef( : , : , 2 ) - Tef( : , : , 3 ) ) )'
% k_s2 * ( 1 - l_s2 / r_sv2 ) * ( r_s2 * ( Tef( : , : , 4 ) - Tef( : , : , 5 ) ) )'
% k_s3 * ( 1 - l_s3 / r_sv3 ) * ( r_s3 * ( Tef( : , : , 6 ) - Tef( : , : , 7 ) ) )'
	
dzt = ( T.' * M * T ) \ ...
    ( T.' * ( fg - M * D * u ) + fj + ...
    NormalForce * ( [ 0 0 1 ] * Tef( : , : , 1 ) )' + ...
    FrictionForce * ( [ 0 -1 0 ] * Tef( : , : , 1 ) )' + ...
    k_s4 * ( 1 - l_s4 / r_sv4 ) * ( r_s4 * ( Tef( : , : , 2 ) - Tef( : , : , 3 ) ) )' + ...
    k_s2 * ( 1 - l_s2 / r_sv2 ) * ( r_s2 * ( Tef( : , : , 4 ) - Tef( : , : , 5 ) ) )' + ...
    k_s3 * ( 1 - l_s3 / r_sv3 ) * ( r_s3 * ( Tef( : , : , 6 ) - Tef( : , : , 7 ) ) )' ...
	) ;

%     ( k_g * r_sgv - c_vg * v_g(3) ) * ( [ 0 0 1 ] * Tef( : , : , 1 ) )' + ...
% NormalForce * meu * sign( v_g(3) ) * ( [ 0 -1 0 ] * Tef( : , : , 1 ) )' + ...
dz = [ u ; dzt ];

