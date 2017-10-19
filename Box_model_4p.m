addpath('~/Dropbox/MOCM/WEILEI/myfunc')
load DB2_May2005.mat
p.dp   = depth;
p.Chl  = Chla;
p.POC  = POC;
p.Phyo = Phyeo;
p.chl  = chla;
p.poc  = poc;
p.phyo = phyeo;
% interpolate nan value in POC
p.Chl_std  = std(p.Chl);
p.chl_std  = std(p.chl);
p.POC_std  = std(p.POC);
p.poc_std  = std(p.poc);
p.Phyo_std = std(p.Phyo);
p.phyo_std = std(p.phyo);
idata = length(p.Chl);
spy      = 365*24*60*60;

depth = p.dp;
jj = length(depth);
M2d = ones(1,jj);

grd = buildgrd(p,M2d);
b  = 0.84;   % Martin Curve exponential
r1 = 1;      % remineralization rate
r2 = 1;      % remineralization rate
p.r3 = 2.1;
a  = 3;      % aggregation rate
d  = 150;    % disaggregation rate
p.eta = 1;     % coefficient to convert conc. to production rate

load Med_xhat_4P.mat
%alpha = linspace(0.1,10,100);
%beta = linspace(20,60,100);
alpha = R.alpha;
beta = R.beta;
[X,Y] = meshgrid(alpha,beta);

logZ = 0*X;
x0 = R.xhat;
				%x0 = [b;r1;r2;r3;a;d];
x0 = log(x0);
nip = length(x0);
for jj = 1:length(X(:))

    p.alpha = X(jj);
    p.beta  = Y(jj);
    L = @(x) neglogpost_4p(x,p,grd,M2d);
    
    options = optimoptions(@fminunc,'Algorithm','trust-region',...
			   'GradObj','on','Hessian','on','Display','off',...
			   'MaxFunEvals',1000,'MaxIter',1000,'TolX',1e-7,...
			   'DerivativeCheck','off','FinDiffType', ...
			   'central','TolFun',1e-7,'PrecondBandWidth',Inf);
    
    [xhat,fval,exitflag] = fminunc(L,x0,options);

    [f,dfdx,d2fdx2] = neglogpost_4p(xhat,p,grd,M2d);
    HH = d2fdx2;
    logZ(jj) = -f-0.5*log(det(HH))+6/2*log(p.alpha)+(6*idata)/2*log(p.beta);
end
 
%figure(2)
%contourf(X,Y,logZ);colorbar
%set(gca,'XTick',[0:1:10])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9','10'})
%set(gca,'YTick',[20:5:60])
%set(gca,'YTickLabel',{'20','25','30','35','40','45','50','55','60'})
%text('Interpreter','latex','String','$$\Lambda = 6.12$$','Position',[35 55])
%text('Interpreter','latex','String','$$\Gamma = 36.97$$','Position',[35 53])
%text(5,55,'MedFlux')
%xlabel('\Lambda-scaling factor for parameter')
%ylabel('\Gamma-scaling factor for data')

imax = find(logZ(:)==max(logZ(:)));
p.alpha = X(imax);
p.beta  = Y(imax);

[xhat,fval,exitflag] = fminunc(L,x0,options);
[f,dfdx,d2fdx2] = neglogpost_4p(xhat,p,grd,M2d);
HH = d2fdx2;

N_p = nip-alpha*(trace(HH)).^(-1);
error = sqrt(diag(inv(HH)));
R.logZ = logZ;
R.upbar = (exp(xhat+error)-exp(xhat));
R.lowbar = (exp(xhat)-exp(xhat-error));
R.xhat = exp(xhat);
R.alpha = p.alpha;
R.beta = p.beta;
%fname = sprintf('Med_xhat_4P');
%save(fname,'R');
