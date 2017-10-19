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

idata = length(p.Chl(2:end));

depth = p.dp;
jj = length(depth);
M2d = ones(1,jj);

grd = buildgrd(p,M2d);
p.eta = 1.2;     % coefficient to convert conc. to production rate

load xhat_log_cons_SV.mat
%alpha = 0;%R.alpha;
%beta  = 1;%R.beta;
alpha = linspace(0.5,1.5,20);
beta = linspace(50,70,20);

[X,Y] = meshgrid(alpha,beta);

logZ = 0*X;
x0 = log([150;R.xhat]); % per day

nip = length(x0);
for jj = 1:length(X(:))

    p.alpha = X(jj);
    p.beta  = Y(jj);
    L = @(x) neglogpost_cons_SV(x,p,grd,M2d);
    
    options = optimoptions(@fminunc,'Algorithm','trust-region',...
			   'GradObj','on','Hessian','on','Display','off',...
			   'MaxFunEvals',1000,'MaxIter',1000,'TolX',1e-12,...
			   'DerivativeCheck','off','FinDiffType', ...
			   'central','TolFun',1e-12,'PrecondBandWidth',Inf);
    
    [xhat,fval,exitflag] = fminunc(L,x0,options);

    [f,dfdx,d2fdx2] = neglogpost_cons_SV(xhat,p,grd,M2d);
    HH = d2fdx2;
    logZ(jj) = -f-0.5*log(det(HH))+nip/2*log(p.alpha)+(6*idata)/2*log(p.beta);
end

%figure(2)
%contourf(X,Y,logZ);colorbar
%set(gca,'XTick',[0.1:0.1:1])
%set(gca,'XTickLabel',{'1','0.10','0.20','0.30','0.40','0.50'})
%set(gca,'YTick',[25:1:45])
%set(gca,'YTickLabel',{'250','275','300','325','350'})
%text('Interpreter','latex','String','$$\Lambda = 0.23$$','Position',[1 230])
%text('Interpreter','latex','String','$$\Gamma = 312.84$$','Position',[1.5 225])
%text(0.2,330,'BarFlux3')
%xlabel('\Lambda-scaling factor for parameter')
%ylabel('\Gamma-scaling factor for data')

imax = find(logZ(:)==max(logZ(:)));
p.alpha = X(imax);
p.beta  = Y(imax);
		   % recalcualte xhat based on optimal alpha and beta;
L = @(x) neglogpost_cons_SV(x,p,grd,M2d);
[xhat,fval,exitflag] = fminunc(L,x0,options);

[f,dfdx,d2fdx2] = neglogpost_cons_SV(xhat,p,grd,M2d);
HH = d2fdx2;
N_p = nip-alpha*trace(inv(HH));

error = sqrt(diag(inv(HH)));
R.upbar = (exp(xhat+error)-exp(xhat));
R.lowbar = (exp(xhat)-exp(xhat-error));
R.xhat = exp(xhat);
R.alpha = p.alpha;
R.beta = p.beta;
R.logZ = logZ;
%fname = sprintf('xhat_log_cons_SV');
%save(fname,'R');
