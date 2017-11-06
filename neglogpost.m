%+++++ function [f,dfdx,d2fdx2] = neglogpost(x,p,grd,M2d)
% this function is used to calculated objective function
% value(f), along with gradient(dfdx) and hessian(d2fdx2) 
% towards parameters. x is a vector of parameters, p is a 
% structure that contains data and constants.M2d is a mask.
function [f,dfdx,d2fdx2] = neglogpost(x,p,grd,M2d)

nip    = length(x);
d2fdx2 = zeros(nip,nip);
dx = sqrt(-1)*eps.^3*eye(nip);
% prescribe prior and variance.
prior = [-0.11;-1.15;-0;-0;1.70;5.85];
U = d0([1/0.16;1/2.98;1/3.98;1/3.98;1/4.00;1/6.64]);
alpha = p.alpha;
beta  = p.beta;

% the gradient is calculated analytically.
% the for loop is used to calculated Hessian matrix by 
% using complex step method. 
for ii = 1:nip
  x = real(x)+dx(:,ii);
  p.b  = exp(x(1));
  p.r1 = exp(x(2));
  p.r2 = exp(x(3));
  p.r3 = exp(x(4));
  p.a  = exp(x(5));
  p.d  = exp(x(6));
  
  [PFD,dPFDdb,dPFDdd] = buildPFD_v2(M2d,p,grd);

  iocn = find(M2d(:));
  tmp = M2d;
  tmp(:,1) = 0;
  ib = find(tmp(iocn));
  
  [M,D] = Pcycle(p,PFD,dPFDdb,dPFDdd,M2d);

  ep = (x-prior);
  W  = speye(length(iocn)-1);
  jj = length(ib);
  
  POC_mod = [M(0*jj+1:1*jj)];
  poc_mod = [M(1*jj+1:2*jj)];
  Chl_mod = [M(2*jj+1:3*jj)];
  chl_mod = [M(3*jj+1:4*jj)];
  Phy_mod = [M(4*jj+1:5*jj)];
  phy_mod = [M(5*jj+1:end)];
  
  e1 = (POC_mod-p.POC(2:end));
  e2 = (poc_mod-p.poc(2:end));
  e3 = (Chl_mod-p.Chl(2:end));
  e4 = (chl_mod-p.chl(2:end));
  e5 = (Phy_mod-p.Phyo(2:end));
  e6 = (phy_mod-p.phyo(2:end));

  f1 = 0.5*(e1.'*W*e1)/p.POC_std;
  f2 = 0.5*(e2.'*W*e2)/p.poc_std;
  f3 = 0.5*(e3.'*W*e3)/p.Chl_std;
  f4 = 0.5*(e4.'*W*e4)/p.chl_std;
  f5 = 0.5*(e5.'*W*e5)/p.Phyo_std;
  f6 = 0.5*(e6.'*W*e6)/p.phyo_std;
  f = f1+f2+f3+f4+f5+f6;
  f = beta*f+0.5*alpha*(ep.'*U*ep);
  %dfdx = e.'*W*dedx;
  % dedx = dedC*dCdx;
  % dCdx = dCdp*dpdx;
  dpdx = diag(exp(x));
  
  de1dx = D(1:jj,:)*dpdx;
  df1dx = (e1.'*W*de1dx)/p.POC_std;

  de2dx = D(jj+1:2*jj,:)*dpdx;
  df2dx = (e2.'*W*de2dx)/p.poc_std;

  de3dx = D(2*jj+1:3*jj,:)*dpdx;
  df3dx = (e3.'*W*de3dx)/p.Chl_std;

  de4dx = D(3*jj+1:4*jj,:)*dpdx;
  df4dx = (e4.'*W*de4dx)/p.chl_std;

  de5dx = D(4*jj+1:5*jj,:)*dpdx;
  df5dx = (e5.'*W*de5dx)/p.Phyo_std;

  de6dx = D(5*jj+1:6*jj,:)*dpdx;
  df6dx = (e6.'*W*de6dx)/p.phyo_std;

  dfdx = df1dx+df2dx+df3dx+df4dx+df5dx+df6dx;
  dfdx = beta*(dfdx)+alpha*ep.'*U;

  dfdx_test(ii) = imag(f)./eps.^3;
  d2fdx2(ii,:) = imag(dfdx)./eps.^3;

end

f = real(f);
dfdx = real(dfdx);

O = [p.POC(2:end);p.poc(2:end);p.Chl(2:end);...
     p.chl(2:end);p.Phyo(2:end);p.phyo(2:end)];

%%%%%%%% Comment this out of model versus observation plot %%%%%%%%%%
%figure(1)
%loglog(real(M),O,'*',[1e-8:10],[1e-8:10])

%xlim([1e-8,10]);
%ylim([1e-8,10]);
%r2 = rsquare(real(M),O);
%txt = sprintf('R^2 = %.2f',r2);
%text(0.01,1,txt)
%xlabel('Model prediction (\mumol L^-^1)')
%ylabel('Observation (\mumol L^-^1)')
%keyboard
%%%%%%%% Comment this out of model versus observation plot %%%%%%%%%%
%++++++++++++++++++++++++++++++++++++++++++++++


%%%%%%++++++function [M,D] = Pcycle(p,PFD,dPFDdb,dPFDdd,M2d)
% function Pcycle output model field (M) and first derivative(D) 
function [M,D] = Pcycle(p,PFD,dPFDdb,dPFDdd,M2d)
  
  r1 = p.r1;
  r2 = p.r2;
  r3 = p.r3;
  a  = p.a;
  d  = p.d;
  eta = p.eta;
  
  jj = length(p.POC);
  I = speye(jj);
  iocn = find(M2d(:));
  tmp = M2d;
  tmp(:,2:end) = 0;
  iu = find(tmp(iocn));
  tmp = M2d;
  tmp(:,1) = 0;
  ib = find(tmp(iocn));
  
  L11 = d*I+PFD;       %dPOC_ldt
  L12 = -a*I;          %dPOC_ldt
  L21 = -d*I;          %dPOC_sdt
  L22 = (a+r1)*I;      %dPOC_sdt
  L33 = d*I+PFD;       %dChl_ldt
  L34 = -a*I;          %dChl_ldt
  L43 = -d*I;          %dChl_sdt
  L44 = (a+r3)*I;      %dChl_sdt
  L55 = d*I+PFD;       %dPhe_ldt
  L56 = -a*I;          %dPhe_ldt
  L64 = -r3*I;         %dPhe_sdt
  L65 = -d*I;          %dPhe_sdt
  L66 = (a+r2)*I;      %dPhe_sdt
  
  LI = 0*I;
  
  LIBB  = LI(ib,:);   LIBB = LIBB(:,ib);
  L11BB = L11(ib,:); L11BB = L11BB(:,ib);
  L12BB = L12(ib,:); L12BB = L12BB(:,ib);
  L21BB = L21(ib,:); L21BB = L21BB(:,ib);
  L22BB = L22(ib,:); L22BB = L22BB(:,ib);
  L44BB = L44(ib,:); L44BB = L44BB(:,ib);
  L64BB = L64(ib,:); L64BB = L64BB(:,ib);
  L66BB = L66(ib,:); L66BB = L66BB(:,ib);
  L33BB = L11BB; L34BB = L12BB; L43BB = L21BB;
  L55BB = L11BB; L56BB = L12BB; L65BB = L21BB;
  
  LIBU  = LI(ib,:);  LIBU  = LIBU(:,iu);
  L11BU = L11(ib,:); L11BU = L11BU(:,iu);
  L12BU = L12(ib,:); L12BU = L12BU(:,iu);
  L21BU = L21(ib,:); L21BU = L21BU(:,iu);
  L22BU = L22(ib,:); L22BU = L22BU(:,iu);
  L44BU = L44(ib,:); L44BU = L44BU(:,iu);
  L64BU = L64(ib,:); L64BU = L64BU(:,iu);
  L66BU = L66(ib,:); L66BU = L66BU(:,iu);
  L33BU = L11BU; L34BU = L12BU; L43BU = L21BU;
  L55BU = L11BU; L56BU = L12BU; L65BU = L21BU;
  
  
  JBB = [[L11BB, L12BB,  LIBB,  LIBB,  LIBB,   LIBB];...%POC_l
	 [L21BB, L22BB,  LIBB,  LIBB,  LIBB,   LIBB];...%POC_s
	 [ LIBB,  LIBB, L33BB, L34BB,  LIBB,   LIBB];...%Chl_l
	 [ LIBB,  LIBB, L43BB, L44BB,  LIBB,   LIBB];...%Chl_s
	 [ LIBB,  LIBB,  LIBB,  LIBB, L55BB,  L56BB];...%Phe_l
	 [ LIBB,  LIBB,  LIBB, L64BB, L65BB,  L66BB]];  %Phe_s
  
  JBU = [[L11BU, L12BU,  LIBU,  LIBU,  LIBU,  LIBU];...
	 [L21BU, L22BU,  LIBU,  LIBU,  LIBU,  LIBU];...
	 [ LIBU,  LIBU, L33BU, L34BU,  LIBU,  LIBU];...
	 [ LIBU,  LIBU, L43BU, L44BU,  LIBU,  LIBU];...
	 [ LIBU,  LIBU,  LIBU,  LIBU, L55BU, L56BU];...
	 [ LIBU,  LIBU,  LIBU, L64BU, L65BU, L66BU]];
  
  poc = M2d;chl = M2d;phyo = M2d;
  POC = M2d;Chl = M2d;Phyo = M2d;
  poc(iocn) = p.poc; chl(iocn) = p.chl; phyo(iocn) = p.phyo; 
  POC(iocn) = p.POC; Chl(iocn) = p.Chl; Phyo(iocn) = p.Phyo;
  
  OBS = [POC(iocn(iu)); poc(iocn(iu)); Chl(iocn(iu));...
	 chl(iocn(iu));Phyo(iocn(iu));phyo(iocn(iu))];
  
  FM     = mfactor(JBB);
  rhs    = -JBU*OBS*eta;
  M      = mfactor(FM,rhs);
  N      = OBS*eta;
  
  dL11dd = I+dPFDdd;    %dPOC_ldt
  dL11db = dPFDdb;      %dPOC_ldt
  dL12da = -I;          %dPOC_ldt
  dL21dd = -I;          %dPOC_sdt
  dL22da =  I;          %dPOC_sdt
  dL22dr1 = I;          %dPOC_sdt
  dL33dd = I+dPFDdd;    %dChl_ldt
  dL33db = dPFDdb;      %dChl_ldt
  dL34da = -I;          %dChl_ldt
  dL43dd = -I;          %dChl_sdt
  dL44da = I;           %dChl_sdt
  dL44dr3 = I;          %dChl_sdt
  dL55dd = I+dPFDdd;    %dPhe_ldt
  dL55db = dPFDdb;      %dPhe_ldt
  dL56da = -I;          %dPhe_ldt
  dL64dr3 = -I;         %dPhe_sdt
  dL65dd = -I;          %dPhe_sdt
  dL66da = I;           %dPhe_sdt
  dL66dr2 = I;          %dPhe_sdt

  % dJdd -----------------
  dL11BBdd = dL11dd(ib,:);
  dL11BBdd = dL11BBdd(:,ib);
  dL33BBdd = dL11BBdd;
  dL55BBdd = dL11BBdd;

  dL21BBdd = dL21dd(ib,:);
  dL21BBdd = dL21BBdd(:,ib);
  dL43BBdd = dL21BBdd;
  dL65BBdd = dL21BBdd;

  dL11BUdd = dL11dd(ib,:);
  dL11BUdd = dL11BUdd(:,iu);  
  dL33BUdd = dL11BUdd;
  dL55BUdd = dL11BUdd;

  dL21BUdd = dL21dd(ib,:);
  dL21BUdd = dL21BUdd(:,iu);
  dL43BUdd = dL21BUdd;
  dL65BUdd = dL21BUdd;
    
  dJBBdd = [[dL11BBdd, LIBB,  LIBB,     LIBB,  LIBB,      LIBB];...%POC_l
	    [dL21BBdd, LIBB,  LIBB,     LIBB,  LIBB,      LIBB];...%POC_s
	    [ LIBB,    LIBB,  dL33BBdd, LIBB,  LIBB,      LIBB];...%Chl_l
	    [ LIBB,    LIBB,  dL43BBdd, LIBB,  LIBB,      LIBB];...%Chl_s
	    [ LIBB,    LIBB,  LIBB,     LIBB,  dL55BBdd,  LIBB];...%Phe_l
	    [ LIBB,    LIBB,  LIBB,     LIBB,  dL65BBdd,  LIBB]];  %Phe_s
  
  dJBUdd = [[dL11BUdd, LIBU,  LIBU,     LIBU,  LIBU,     LIBU];...
	    [dL21BUdd, LIBU,  LIBU,     LIBU,  LIBU,     LIBU];...
	    [    LIBU, LIBU,  dL33BUdd, LIBU,  LIBU,     LIBU];...
	    [    LIBU, LIBU,  dL43BUdd, LIBU,  LIBU,     LIBU];...
	    [    LIBU, LIBU,  LIBU,     LIBU,  dL55BUdd, LIBU];...
	    [    LIBU, LIBU,  LIBU,     LIBU,  dL65BUdd, LIBU]];


  % dJdb-------------------
  dL11BBdb = dL11db(ib,:); dL11BBdb = dL11BBdb(:,ib);
  dL33BBdb = dL11BBdb;
  dL55BBdb = dL11BBdb;

  dL11BUdb = dL11db(ib,:); dL11BUdb = dL11BUdb(:,iu);  
  dL33BUdb = dL11BUdb;
  dL55BUdb = dL11BUdb;
    
  dJBBdb = [[dL11BBdb, LIBB,  LIBB,     LIBB,  LIBB,      LIBB];...%POC_l
	    [ LIBB,    LIBB,  LIBB,     LIBB,  LIBB,      LIBB];...%POC_s
	    [ LIBB,    LIBB,  dL33BBdb, LIBB,  LIBB,      LIBB];...%Chl_l
	    [ LIBB,    LIBB,  LIBB,     LIBB,  LIBB,      LIBB];...%Chl_s
	    [ LIBB,    LIBB,  LIBB,     LIBB,  dL55BBdb,  LIBB];...%Phe_l
	    [ LIBB,    LIBB,  LIBB,     LIBB,  LIBB,      LIBB]];  %Phe_s
  
  dJBUdb = [[dL11BUdb, LIBU,  LIBU,     LIBU,  LIBU,     LIBU];...
	    [    LIBU, LIBU,  LIBU,     LIBU,  LIBU,     LIBU];...
	    [    LIBU, LIBU,  dL33BUdb, LIBU,  LIBU,     LIBU];...
	    [    LIBU, LIBU,  LIBU,     LIBU,  LIBU,     LIBU];...
	    [    LIBU, LIBU,  LIBU,     LIBU,  dL55BUdb, LIBU];...
	    [    LIBU, LIBU,  LIBU,     LIBU,  LIBU,     LIBU]];

  % dJda ---------------------
  dL12BBda = dL12da(ib,:); dL12BBda = dL12BBda(:,ib);
  dL34BBda = dL12BBda;
  dL56BBda = dL12BBda;
  dL22BBda = dL22da(ib,:); dL22BBda = dL22BBda(:,ib);  
  dL44BBda = dL22BBda;
  dL66BBda = dL22BBda;

  dL12BUda = dL12da(ib,:); dL12BUda = dL12BUda(:,iu);
  dL34BUda = dL12BUda;
  dL56BUda = dL12BUda;
  dL22BUda = dL22da(ib,:); dL22BUda = dL22BUda(:,iu);  
  dL44BUda = dL22BUda;
  dL66BUda = dL22BUda;

  dJBBda = [[LIBB, dL12BBda,  LIBB,     LIBB,  LIBB,      LIBB];...%POC_l
	    [LIBB, dL22BBda,  LIBB,     LIBB,  LIBB,      LIBB];...%POC_s
	    [LIBB,     LIBB,  LIBB, dL34BBda,  LIBB,      LIBB];...%Chl_l
	    [LIBB,     LIBB,  LIBB, dL44BBda,  LIBB,      LIBB];...%Chl_s
	    [LIBB,     LIBB,  LIBB,     LIBB,  LIBB,  dL56BBda];...
	    [LIBB,     LIBB,  LIBB,     LIBB,  LIBB,  dL66BBda]];  %Phe_s
  
  dJBUda = [[LIBU, dL12BUda,  LIBU,     LIBU,  LIBU,     LIBU];...
	    [LIBU, dL22BUda,  LIBU,     LIBU,  LIBU,     LIBU];...
	    [LIBU,     LIBU,  LIBU, dL34BUda,  LIBU,     LIBU];...
	    [LIBU,     LIBU,  LIBU, dL44BUda,  LIBU,     LIBU];...
	    [LIBU,     LIBU,  LIBU,     LIBU,  LIBU, dL56BUda];...
	    [LIBU,     LIBU,  LIBU,     LIBU,  LIBU, dL66BUda]];
  
  % dJdr1 ---------------
  dL22BBdr1 = dL22dr1(ib,:); dL22BBdr1 = dL22BBdr1(:,ib);
  dL22BUdr1 = dL22dr1(ib,:); dL22BUdr1 = dL22BUdr1(:,iu);
  
  dJBBdr1 = [[LIBB,      LIBB,  LIBB, LIBB, LIBB, LIBB];...%POC_l
	     [LIBB, dL22BBdr1,  LIBB, LIBB, LIBB, LIBB];...%POC_s
	     [LIBB,      LIBB,  LIBB, LIBB, LIBB, LIBB];...%Chl_l
	     [LIBB,      LIBB,  LIBB, LIBB, LIBB, LIBB];...%Chl_s
	     [LIBB,      LIBB,  LIBB, LIBB, LIBB, LIBB];...%Phe_l
	     [LIBB,      LIBB,  LIBB, LIBB, LIBB, LIBB]];  %Phe_s
  
  dJBUdr1 = [[LIBU,      LIBU, LIBU, LIBU, LIBU, LIBU];...
	     [LIBU, dL22BUdr1, LIBU, LIBU, LIBU, LIBU];...
	     [LIBU,      LIBU, LIBU, LIBU, LIBU, LIBU];...
	     [LIBU,      LIBU, LIBU, LIBU, LIBU, LIBU];...
	     [LIBU,      LIBU, LIBU, LIBU, LIBU, LIBU];...
	     [LIBU,      LIBU, LIBU, LIBU, LIBU, LIBU]];


  % dJdr2 ----------------
  dL66BBdr2 = dL66dr2(ib,:); dL66BBdr2 = dL66BBdr2(:,ib);
  dL66BUdr2 = dL66dr2(ib,:); dL66BUdr2 = dL66BUdr2(:,iu);
  
  dJBBdr2 = [[LIBB, LIBB, LIBB, LIBB,  LIBB,      LIBB];...%POC_l
	     [LIBB, LIBB, LIBB, LIBB,  LIBB,      LIBB];...%POC_s
	     [LIBB, LIBB, LIBB, LIBB,  LIBB,      LIBB];...%Chl_l
	     [LIBB, LIBB, LIBB, LIBB,  LIBB,      LIBB];...%Chl_s
	     [LIBB, LIBB, LIBB, LIBB,  LIBB,      LIBB];...%Phe_l
	     [LIBB, LIBB, LIBB, LIBB,  LIBB, dL66BBdr2]];  %Phe_s
  
  dJBUdr2 = [[LIBU, LIBU, LIBU, LIBU,  LIBU,      LIBU];...
	     [LIBU, LIBU, LIBU, LIBU,  LIBU,      LIBU];...
	     [LIBU, LIBU, LIBU, LIBU,  LIBU,      LIBU];...
	     [LIBU, LIBU, LIBU, LIBU,  LIBU,      LIBU];...
	     [LIBU, LIBU, LIBU, LIBU,  LIBU,      LIBU];...
	     [LIBU, LIBU, LIBU, LIBU,  LIBU, dL66BUdr2]];


  % dJdr3 ---------------
  dL44BBdr3 = dL44dr3(ib,:); dL44BBdr3 = dL44BBdr3(:,ib);
  dL44BUdr3 = dL44dr3(ib,:); dL44BUdr3 = dL44BUdr3(:,iu);

  dL64BBdr3 = dL64dr3(ib,:); dL64BBdr3 = dL64BBdr3(:,ib);
  dL64BUdr3 = dL64dr3(ib,:); dL64BUdr3 = dL64BUdr3(:,iu);
  
  dJBBdr3 = [[LIBB, LIBB,  LIBB, LIBB,      LIBB, LIBB];...%POC_l
	     [LIBB, LIBB,  LIBB, LIBB,      LIBB, LIBB];...%POC_s
	     [LIBB, LIBB,  LIBB, LIBB,      LIBB, LIBB];...%Chl_l
	     [LIBB, LIBB,  LIBB, dL44BBdr3, LIBB, LIBB];...%Chl_s
	     [LIBB, LIBB,  LIBB, LIBB,      LIBB, LIBB];...%Phe_l
	     [LIBB, LIBB,  LIBB, dL64BBdr3, LIBB, LIBB]];  %Phe_s
  
  dJBUdr3 = [[LIBU, LIBU,  LIBU,      LIBU,  LIBU, LIBU];...
	     [LIBU, LIBU,  LIBU,      LIBU,  LIBU, LIBU];...
	     [LIBU, LIBU,  LIBU,      LIBU,  LIBU, LIBU];...
	     [LIBU, LIBU,  LIBU, dL44BUdr3,  LIBU, LIBU];...
	     [LIBU, LIBU,  LIBU,      LIBU,  LIBU, LIBU];...
	     [LIBU, LIBU,  LIBU, dL64BUdr3   LIBU, LIBU]];
  
  
  D = -mfactor(FM,[dJBUdb*N+dJBBdb*M,...
		  dJBUdr1*N+dJBBdr1*M,...
		  dJBUdr2*N+dJBBdr2*M,...
		  dJBUdr3*N+dJBBdr3*M,...
		  dJBUda*N+dJBBda*M,...
		  dJBUdd*N+dJBBdd*M]);
