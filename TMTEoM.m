%% TMT EOM Derivator:
% ===========================================
% [ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn ,  Tef , Tcn , Dcn , qf , uf ] = ...
%    TMTEoM ( lc , m , I , j , jkd , g , quatSlct , mtdSlct );
%
% ===========================================
% Author:
%   S.M.Hadi Sadati
%   PhD student - King's College London
%   2014
%
% ===========================================
% Help:
%
% Inputs:
% ( n ) Number of links (every thing that has mass and/or rotational inertia
% ( nq ) Number of generalized coordinates
% lc : nx5
%     first 3 elements are vector of link COM (4th: 0) / external force (4th: 1) / constraint (4th: 2) positions in link frame (1x3),
%     	4th element is indicator for types
%     	and 5th element is the number of previous joint ( 0 for serial mechanism ).
% 	  Deal with links first, then forces and then constraints.
%     Joint (link) frame origin attached to joint
% m : nx1 vector of link masses
% I : 3x3xn cube of link inertia matrices w.r.t. link frame
% j : rx5xn cube of each joint (translation and rotation leading to each joint)
%     translation from previous joint position:
%       ( r , 3:5 , n ) = linear transformation value along previous link frame x-y-z axes 
%     and rotation of new link frame w.r.t. previous link frame
%       ( r , 1 , n ) = 1 : 3 axis indicator number for Euler post-multiplication rotations about sequencive temporary frames' x-y-z axes ),
%     and the rotation values:
%       ( r , 2 , n ) = value.
%     PLACE "inf" FOR "VALUE" IF IT IS ONE OF THE GENERALIZED COORDINATES!!!
% jkd : 3x2xnq cube generalized coordinates 3 pair of 
%     [ spring coeff.s ; viscous damping coeff.s ; external input ]
%     and their initial pos.s ( Only spring coeff. needs it! )
% g : 1x3 gravity acceleration vector
% quatSlct : Select quaternion (1) or matrix transformation (2)
% mtdSlct : Select TMT method (1) or Lagrange Method (2)
%
% Outputs:
% M : Mass matrix
% T : Transformation matrix
% Dd : Damping/Stiffness matrix
% fg : Gravity force virtual work
% fj : Joint stiffness/damping virtual work acting directly on generalized coordinates
% rj : 3 first elements are for Joint absolute position vector in base frame
%	  and the last element the number of the previous element ( 0 for series link part )
% rc : Joint COM absolute position vector in base frame
% vc : Joint COM absolute linear velocity vector in base frame
% wc : Joint COM absolute rotational velocity vector in link frame
% ref: External force/torque position vector
% rcn: Constraint location ( Useful for parallel mechanisms )
% Tef: External force/torque position vector Jacobian w.r.t. generalized coord.s
% Tcn: Constraint location Jacobian w.r.t. generalized coord.s
% Dcn: Remain terms in Constraint 2nd time location derivation
% qf : Generalized coordinates
% uf : Generalized coordinates derivatives
% 
% ============================================
% Example:
%
% syms lc1 lc2 lc3 m1 m2 m3 I1 I2 I3 l1 l2 k1 k2 k3 gx gy gz
% lc = [ lc1 0 0 0 0 ; lc2 0 0 0 0 ; lc3 0 0 0 0 ];
% m = [ m1 , m2 , m3 ];
% I = sym ( zeros ( 3 , 3 , 3 ) );
% I(:,:,1) = I1 * eye ( 3 ); I(:,:,2) = I2 * eye ( 3 ); I(:,:,3) = I3 * eye ( 3 );
% j = sym ( zeros ( 1 , 5 , 3 ) );
% j(:,:,1) = [ 2 inf 0 0 0 ]; j(:,:,2) = [ 2 inf l1 0 0 ]; j(:,:,3) = [ 2 inf l2 0 0 ];
% jkd = sym (zeros ( 3 , 2 , 3 ) );
% jkd(1,:,1) = [ k1 0 ]; jkd(1,:,2) = [ k2 0 ]; jkd(1,:,3) = [ k3 0 ];
% g = [ gx , gy , gz ];
%
% ============================================



%% Main:
function [ M , T , Dd , fg , fj , rj , rc , vc , wc , ref , rcn , Tef , Tcn , Dcn , qf , uf ] = ...
    TMTEoM ( lc , m , I , j , jkd , g , quatSlct , mtdSlct )


% >> Initialization:

lc = sym ( lc );
m = sym ( m );
I = sym ( I );
j = sym ( j );
jkd = sym ( jkd );
g = sym ( g );

[ n , tmp ] = size ( lc ); % number of links and ext. forces and const. positions
nl = 0; % number of the links
nef = 0; % number of ext. forces/torques
ncn = 0; % number of constraints

for i = 1 : n

	switch lc(i,4)
	
		case 0
			nl = nl + 1;
			
		case 1
			nef = nef + 1;
			
		case 2
			ncn = ncn + 1;
			
	end
	
end

[ rl , tmp , r2 ] = size ( j ); % maximum number of rotation in each joint
nq = 0; % number of generalized coordinates
for i = 1 : r2
    for i1 = 1 : rl
        for i2 = 2 : 5
			if j( i1 , i2 , i ) == inf && jkd(4,1,nq+1) == 0
				nq = nq + 1;
			end; end; end; end

M = sym ( zeros ( 6 * nl , 6 * nl ) );
T = sym ( zeros ( 6 * nl , nq ) );
Dd = sym ( zeros ( 6 * nl , nq ) );
fg = sym ( zeros ( 6 * nl , 1 ) );
fj = sym ( zeros ( nq , 1 ) );
rj = sym ( zeros ( nl , 4 ) );
rc = sym ( zeros ( nl , 3 ) );
vc = sym ( zeros ( nl , 3 ) );
wcb = sym ( zeros ( nl , 3 ) );
wc = sym ( zeros ( nl , 3 ) );
Wc = sym ( zeros ( 3 , 3 , nl ) );
ref = sym ( zeros ( nef , 3 ) );
rcn = sym ( zeros ( ncn , 3 ) );
Tef = sym ( zeros ( 3 , nq , nef ) );
Tcn = sym ( zeros ( 3 , nq , ncn ) );
Dcn = sym ( zeros ( 3 , nq , ncn ) );


syms q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16 q17 q18 q19 q20 q21
syms u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21
q = [ q1 q2 q3 q4 q5 q6 q7 q8 q9 q10 q11 q12 q13 q14 q15 q16 q17 q18 q19 q20 q21 ];
u = [ u1 u2 u3 u4 u5 u6 u7 u8 u9 u10 u11 u12 u13 u14 u15 u16 u17 u18 u19 u20 u21 ];
iq = 0;
qf = sym ( [] );
uf = sym ( [] );

R = sym ( zeros ( 3 , 3 , n ) );
TR = sym ( zeros ( 4 , 4 , n ) );
Rt = sym ( zeros ( 3 , 3 , n ) );
TRt = sym ( zeros ( 4 , 4 , n ) );
fgt = sym ( [] );


% >> Derivation:

il = 0; % link counter
ief = 0; % external force counter
icn = 0; % constraint position counter

for i = 1 : n
    
	if lc(i,4) == 0
        
        il = il + 1;
	
		im = 3 * il;
		iI = 3 * nl + 3 * il;
    
		M ( im-2 : im , im-2 : im ) = m(il) * eye ( 3 );
		M ( iI-2 : iI , iI-2 : iI ) = I(:,:,il);
    
	end
	
    R(:,:,i) = sym ( eye ( 3 ) );
    TR(:,:,i) = sym ( eye ( 4 ) );
    
	if lc(i,5) == 0
		for i1 = 1 : rl
        
		
			if j(i1,:,i) == 0
				break; end
        
			for i2 = 2 : 5            
				if j(i1,i2,i) == inf
                
					if jkd(4,1,iq+1) == 0                
						iq = iq + 1;
						iqc = iq;
						qf = [ qf , q(iq) ];
						uf = [ uf , u(iq) ];
					else
						iqc = jkd(4,1,iq+1);
					end
                
					j(i1,i2,i) = q(iqc);
                
					fj(iqc,1) = fj(iqc,1) + ...
						- jkd(1,1,iq) * ( q(iqc) - jkd(1,2,iq) ) ... % spring
						- jkd(2,1,iq) * u(iqc) ... % viscous
						+ jkd(3,1,iq); % external input
                                
				end
			end
        
			[ tmp1 , tmp2 ] = TRm ( j(i1,1:5,i) );
			R(:,:,i) = R(:,:,i) * tmp1;
			TR(:,:,i) = TR(:,:,i) * tmp2;
        
                
		end
	end
    
	
    if i > 1
	
		if lc(i,5) == 0
		
			Rt(:,:,i) = Rt(:,:,i-1) *  R(:,:,i);
			TRt(:,:,i) = TRt(:,:,i-1) *  TR(:,:,i);
			
		else
		
			Rt(:,:,i) = Rt(:,:,lc(i,5)) *  R(:,:,i);
			TRt(:,:,i) = TRt(:,:,lc(i,5)) *  TR(:,:,i);
			
		end
		
		
    else
	
        Rt(:,:,i) = R(:,:,i);
        TRt(:,:,i) = TR(:,:,i);
		
		
    end
	
    
	switch lc(i,4)
	
		case 0
			
			tmp3 = sym ( [] );
			for i1 = 1 : 3
				tmp4 = jacobian ( Rt(:,i1,i) , qf ) * uf.';
				tmp3 = [ tmp3 , tmp4 ];
			end
			Wc(:,:,il) = tmp3 * Rt(:,:,i).'; % angular velocity tensor
			wcb(il,:) = [ Wc(3,2,il) , Wc(1,3,il) , Wc(2,1,il) ]; % w in base frame
			wc(il,:) = ( Rt(:,:,il).' * wcb(il,:).' ).'; % w in link frame
    
            
			for i1 = 1 : 3
				for i2 = 1 : iq
					[ tmp5 , tmp6 ] = coeffs ( wc(il,i1) , uf(i2) );
					if isempty ( tmp5 )
						tmp5 = sym ( 0 );
					end
					T ( iI - 3 + i1 , i2 ) = tmp5(1);
				end; end
    
			Dd ( iI - 2 : iI , 1 : iq ) = jacobian ( wc(il,:).' , qf );
		
			rj(il,1:3) = TRt(1:3,4,il).';
			rj(il,4) = lc(i,5);
			tmp7 = TRt(:,:,i) * [ lc(i,1:3) , 1 ].';
			rc(il,:) = tmp7(1:3).';
			tmp8 = jacobian ( tmp7(1:3) , qf );
			vc(il,:) = ( tmp8 * uf.' ).';
    
			T ( im - 2 : im , 1 : iq ) = tmp8;
			Dd ( im - 2 : im , 1 : iq ) = jacobian ( vc(il,:).' , qf );
    
			fgt = [ fgt , g ];
			
			
		case 1
		
			ief = ief + 1;
			
			tmp9 = TRt(:,:,i) * [ lc(i,1:3) , 1 ].';
			ref(ief,:) = tmp9(1:3).';
			Tef ( : , 1 : iq , ief ) = jacobian ( tmp9(1:3) , qf );
		
		
		case 2
		
			icn = icn + 1;
			
			tmp10 = TRt(:,:,i) * [ lc(i,1:3) , 1 ].';
			rcn(icn,:) = tmp10(1:3).';
			tmp11 = jacobian ( tmp10(1:3) , qf );
			vcn = ( tmp11 * uf.' ).';
    
			Tcn ( : , 1 : iq , icn ) = tmp11;
			Dcn ( : , 1 : iq , icn ) = jacobian ( vcn.' , qf );
		
	end
    
end


for i1 = 1 : n
    
    if lc(i1,5) ~= 0
            
        if lc( i1 - 1 ,4) == 0
            
            il = il + 1;
            rj(il,1:3) = rc( i1 - 1 ,:);
            rj(il,4) = i1 - 1;
            
        else
            
            break
            
        end
        
    end
    
end

fg ( 1 : 3 * nl , 1 ) = M ( 1 : 3 * nl , 1 : 3 * nl ) * fgt.';

if nef == 0; ref = sym ( 0 ); Tef = sym ( 0 ); end % for system without const. or ext. forces
if ncn == 0; rcn = sym ( 0 ); Tcn = sym ( 0 ); Dcn = sym ( 0 ); end


% value = {0*M; 0*T; 0*Dd; 0*fg; 0*fj; 0*rj; 0*rc; 0*vc; 0*wc; 0*ref; 0*rcn; 0*Tef; 0*Tcn; 0*Dcn};
% out = struct('f',value);
% 
% % matlabpool('open', 4);
% parfor i = 1 : 14
%     
%     switch i
%         
%         case 1
% %             out(i).f = simple ( M );
% %             ccode ( M , 'file' , 'M.cpp' ); % save as c file
%             matlabFunction ( M , 'file' , 'M.m' ); % save as matlab function
% %             M_c = ccode ( M );  % save as c format text file
% %             M_o = fopen ( 'M.cpp' , 'w' );
% %             fwrite ( M_o , M_c );
% %             fclose ( M_o );
%             
%         case 2
% %             out(i).f = simple ( T );
% %             ccode ( T , 'file' , 'T.cpp' );
%             matlabFunction ( T , 'file' , 'T.m' );
% %             T_c = ccode ( T );
% %             T_o = fopen ( 'T.cpp' , 'w' );
% %             fwrite ( T_o , T_c );
% %             fclose ( T_o );
%             
%         case 3
% %             out(i).f = simple ( Dd );
% %             ccode ( Dd , 'file' , 'Dd.cpp' );
%             matlabFunction ( Dd , 'file' , 'Dd.m' );
% %             Dd_c = ccode ( Dd );
% %             Dd_o = fopen ( 'Dd.cpp' , 'w' );
% %             fwrite ( Dd_o , Dd_c );
% %             fclose ( Dd_o );
%             
%         case 4
% %             out(i).f = simple ( fg );
% %             ccode ( fg , 'file' , 'fg.cpp' );
%             matlabFunction ( fg , 'file' , 'fg.m' );
% %             fg_c = ccode ( fg );
% %             fg_o = fopen ( 'fg.cpp' , 'w' );
% %             fwrite ( fg_o , fg_c );
% %             fclose ( fg_o );
%             
%         case 5
% %             out(i).f = simple ( fj );
% %             ccode ( fj , 'file' , 'fj.cpp' );
%             matlabFunction ( fj , 'file' , 'fj.m' );
% %             fj_c = ccode ( fj );
% %             fj_o = fopen ( 'fj.cpp' , 'w' );
% %             fwrite ( fj_o , fj_c );
% %             fclose ( fj_o );
%             
%         case 6
% %             out(i).f = simple ( rj );
% %             ccode ( rj , 'file' , 'rj.cpp' );
%             matlabFunction ( rj , 'file' , 'rj.m' );
% %             rj_c = ccode ( rj );
% %             rj_o = fopen ( 'rj.cpp' , 'w' );
% %             fwrite ( rj_o , rj_c );
% %             fclose ( rj_o );
%             
%         case 7
% %             out(i).f = simple ( rc );
% %             ccode ( rc , 'file' , 'rc.cpp' );
%             matlabFunction ( rc , 'file' , 'rc.m' );
% %             rc_c = ccode ( rc );
% %             rc_o = fopen ( 'rc.cpp' , 'w' );
% %             fwrite ( rc_o , rc_c );
% %             fclose ( rc_o );
%             
%         case 8
% %             out(i).f = simple ( vc );
% %             ccode ( vc , 'file' , 'vc.cpp' );
%             matlabFunction ( vc , 'file' , 'vc.m' );
% %             vc_c = ccode ( vc );
% %             vc_o = fopen ( 'vc.cpp' , 'w' );
% %             fwrite ( vc_o , vc_c );
% %             fclose ( vc_o );
%             
%         case 9
% %             out(i).f = simple ( wc );
% %             ccode ( wc , 'file' , 'wc.cpp' );
%             matlabFunction ( wc , 'file' , 'wc.m' );
% %             wc_c = ccode ( wc );
% %             wc_o = fopen ( 'wc.cpp' , 'w' );
% %             fwrite ( wc_o , wc_c );
% %             fclose ( wc_o );
%             
%         case 10
% %             out(i).f = simple ( ref );
% %             ccode ( ref , 'file' , 'rcn.cpp' );
%             matlabFunction ( ref , 'file' , 'ref.m' );
% %             ref_c = ccode ( ref );
% %             ref_o = fopen ( 'ref.cpp' , 'w' );
% %             fwrite ( ref_o , ref_c );
% %             fclose ( ref_o );
%             
%         case 11
% %             out(i).f = simple ( rcn );
% %             ccode ( rcn , 'file' , 'ref.cpp' );
%             matlabFunction ( rcn , 'file' , 'rcn.m' );
% %             rcn_c = ccode ( rcn );
% %             rcn_o = fopen ( 'rcn.cpp' , 'w' );
% %             fwrite ( rcn_o , rcn_c );
% %             fclose ( rcn_o );
%             
%         case 12
% %             out(i).f = simple ( Tef );
% %             ccode ( Tef , 'file' , 'Tef.cpp' );
%             matlabFunction ( Tef , 'file' , 'Tef.m' );
% %             Tef_c = ccode ( Tef );
% %             Tef_o = fopen ( 'Tef.cpp' , 'w' );
% %             fwrite ( Tef_o , Tef_c );
% %             fclose ( Tef_o );
%             
%         case 13
% %             out(i).f = simple ( Tcn );
% %             ccode ( Tcn , 'file' , 'Tcn.cpp' );
%             matlabFunction ( Tcn , 'file' , 'Tcn.m' );
% %             Tcn_c = ccode ( Tcn );
% %             Tcn_o = fopen ( 'Tcn.cpp' , 'w' );
% %             fwrite ( Tcn_o , Tcn_c );
% %             fclose ( Tcn_o );
%             
%         case 14
% %             out(i).f = simple ( Dcn );
% %             ccode ( Dcn , 'file' , 'Dcn.cpp' );
%             matlabFunction ( Dcn , 'file' , 'Dcn.m' );
% %             Dcn_c = ccode ( Dcn );
% %             Dcn_o = fopen ( 'Dcn.cpp' , 'w' );
% %             fwrite ( Dcn_o , Dcn_c );
% %             fclose ( Dcn_o );
% 
%     end
%     
% end
% 
% M = out(1).f;
% T = out(2).f;
% Dd = out(3).f;
% fg = out(4).f;
% fj = out(5).f;
% rj = out(6).f;
% rc = out(7).f;
% vc = out(8).f;
% wc = out(9).f;
% ref = out(10).f;
% rcn = out(11).f;
% Tef = out(12).f;
% Tcn = out(13).f;
% Dcn = out(14).f;


%% Complementary Functions:

function [ R , TR ] = TRm ( r ) % base rotation matices

i = r(1);
x = r(2);

switch i
    
    case 1
        R = [ ...
            1       0        0 ;
            0  cos(x)  -sin(x);
            0  sin(x)   cos(x)];
        
    case 2
        R = [ ...
             cos(x)  0  sin(x);
                  0  1       0;
            -sin(x)  0  cos(x)];
        
    case 3
        R = [ ...
            cos(x)  -sin(x)  0;
            sin(x)   cos(x)  0;
                 0        0  1];
    otherwise
        R = sym ( eye ( 3 ) );
        
end

TR = sym ( [ [ R , r(3:5).' ] ; 0 0 0 1 ] );
 

%% Notes:
% Each TR contains translation and then a rotation

