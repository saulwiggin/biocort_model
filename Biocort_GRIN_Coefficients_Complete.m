%RAY TRACING PROGRAM TO CALCULATE THE CHROMATIC SEIDEL ABERRATIONS IN A 
%RADIAL GRIN LENS USING QUASI-INVARIANT APPROACH + RAY TRACING. SEE biocort
%LIMITATIONS: RADIAL INDEX ONLY ie.WOOD LENS

% Refractive idnex gradient n(r) = n_0(1-kr.^2)
n = linspace(3.5,5.5,100);
k = 2.4;
g = k.^(-0.5);
d = 0.1;
nsurface = 99;
nrays = 1;

% Surface curvature
a4 = 0.00345;
c = linspace(1,30,1);

% Define abbe number
v_g = 50;
v_h = 50;
dell_lam_k=50;
dell_lam=50;

%Define color abberation 
lambda_lamp = 0;

% Define cheif and marginal ray vectors
u = linspace(1,30,100);
h = linspace(1,1,1);
w = linspace(1,30,100);
m = linspace(1,1,1);

% Define seidel coefficients
Wsph = 0;
Wcoma = 0;
Wastig = 0;
Wcurv = 0;
Wdist = 0;

% Unknown variables:
N4 = 1;

% Refraction at a single surface
for i = 1: nsurface
  
n0 = n(i);
n0_dash = n(i+1);
  
    % Paraxial quasi-invariants
e1 = k.*h.^2+u.^2;
e2 = k.*h.*m+u.*w;
e3 = k.*m.^2+w.^2;
n0i = n0.*h*a4-n0.*u;
n0j = n0.*m*a4-n0.*w;

% Lagrange invariant
H = n.*(u.*m-w.*h);
A = n.*u;
A_bar = n.*w;

    % RAY TRACING EQUATIONS FROM PREVIOUS SURFACE (u, h, w, m)
    % k must be positive. if k<0 see eq 9 biocort. if k=0 standard ray
    % trace equations apply.
    % for the marginal ray
    u_dash = u*cos(g*d)+h*g*sin(g*d);
    h_dash = -(u./g)*sin(g*d)+h*cos(g*d);
    % for the chief ray
    w_dash = w*cos(g*d)+m*g*sin(g*d);
    m_dash = -(w./g)*sin(g*d)+m*cos(g*d);
    
    % calculate abbe number contribution
    dell_lam_dash = (n0_dash-1)./(n0*v_h);
    dell_lam_k_dash = 1./v_g-(n0-1)./(n0*v_h);
    
%setup dell values from previous surface 
dell_h_u = h_dash.*u_dash - h.*u;
dell_h_w = h_dash.*w_dash-h.*w;
dell_lam_n0 = (dell_lam_dash - dell_lam);
dell_lam_k = dell_lam_k_dash - dell_lam_k;
%and for seidel
dell_u_n0 = u_dash./n0_dash-u./n0;
dell_n = n0_dash-n0;
dell_n_k = n0_dash.*k-n0*k;
dell_h_u = h_dash.*u_dash-h.*u;
dell_m_u = m_dash.*u_dash-h.*u;
dell_u_w = u_dash.*w_dash - u.*w;
dell_h_w = h_dash.*w_dash - h.*w;
dell_w2 = w_dash.*w_dash - w.*w;

% surface contibutions to colour
S_lam1 = h.*n0i.*dell_lam_n0;
S_lam2 = h.*n0j.*dell_lam_n0;

%surface contributions to seidel coefficients
%homogeneous surface contributions
Ssph = -(n0.*u).^2.*h.*dell_u_n0+a4.*dell_n.*h.^4;
Scoma = (n0.*u).*(n0.*w).*h.*dell_u_n0;
Sastig = -(n0.*w).^2.*h.*dell_u_n0;
Scurv = -((n0.*w).^2.*h.*dell_u_n0+H.^2.*c.*dell_n);
Sdist = -(n0.*w).^3.*h.*dell_n.^2 + m.*(n0.*w.*(2.*h.*(n.*w)-m.*(n0.*u)).*c.*dell_n);

%inhomogeneous surface contributions 
Ssph_s = -2*h.^4.*c.*dell_n_k;
Scoma_s = -2.*h.^3.*m.*c.*dell_n_k;
Sastig_s = -2.*h.^2.*m.^2.*c.*dell_n_k;
Sdist_s = -2.*h.*m.^3.*c.*dell_n_k;

% equations for lateral and transverse colour

T_lam1 = (n0_dash./2).*(dell_lam_k)*(d*e1+dell_h_u);
T_lam2 = (n0_dash./2).*(dell_lam_k)*(d*e2+dell_h_w);

%and for seidel coefficients
Tsph = n0.*d.*e1.^2.*(1-(3*N4/2))-n*(1+N4).*dell_h_u.^3+5/2*n0*N4*e1.*dell_h_u;
Tcoma= n0.*d.*e1.*e2.*(1-(3.*N4/2))-n0.*(1+N4).*dell_h_u.^2.*w + (5/2).*n0.*N4.*e3.*dell_h_u-2.*N4.*H.*dell_h_u.^2.*w;
Tcurv = k*d*H.^2/n0;
Tastig = n0*d*e2.^2*(1-(3*N4/2))-n0*(1+N4).*dell_h_u.*w.^2 + (5/2)*n0*N4*e3.*dell_h_u-2*N4*H.*dell_u_w-1/2*N4.*Tcurv;
Tdist = n0.*d.*e2.*e3.*(1-(3.*N4/2))-n0.*(1+N4).*dell_h_w.^3+(5/2).*n0.*N4.*e3.*dell_h_w-(1/2).*N4.*H.*dell_w2;

% TOTAL CHROMATIC PARAXIAL ABERRATIONS COEFFICIENTS
% for each surface accumulating each time
C = lambda_lamp + S_lam1 + S_lam2 +T_lam1 + T_lam2

%and for seidel aberrations
Wsph = Wsph+(1/8)*(Tsph + Ssph + Ssph_s)
Wcoma = Wcoma+(1/2)*( Tcoma + Scoma_s + Scoma)
Wastig = Wastig+(1/2)*(Tastig + Sastig_s + Sastig)
Wcurv = Wcurv+(1/4)*(Tastig +Sastig_s + Sastig + Tcurv + Scurv)
Wdist = Wdist+(1/2)*(Tdist + Sdist_s + Sdist)

% %sum all wavefront aberrations for all surfaces
% 
% Wsph_sum = sum(Wsph);
% Wcoma_sum = sum(Wcoma);
% Wastig_sum = sum(Wastig);
% Wcurv_sum = sum(Wcurv);
% Wdist_sum = sum(Wdist)

%store surface values
h = h_dash;
w = w_dash;
u = u_dash;
m = m_dash;
dell_lam = dell_lam_dash;
dell_lam_l = dell_lam_k_dash;
dell_lam_k = dell_lam_k_dash;

end

save results_aberrations Wsph Wcoma Wastig Wcurv Wdist C
