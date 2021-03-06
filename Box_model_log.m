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
idata      = length(p.Chl);

spd        = 365;

depth = p.dp;
jj  = length(depth);
M2d = ones(1,jj);
grd = buildgrd(p,M2d);

p.b   = 0.84;   % Martin Curve exponential
p.r1  = 1;      % remineralization rate
p.r2  = 1;      % remineralization rate
p.r3  = 1;
p.a   = 3;      % aggregation rate
p.d   = 150;    % disaggregation rate
p.eta = 1;     % coefficient to convert conc. to production rate

alpha = linspace(0.1,1,20);
beta  = linspace(0.1,1,20);
%alpha = 0;%R.alpha;
%beta  = 1;%R.beta;
[X,Y] = meshgrid(alpha,beta);

logZ = 0*X;

x0 = [p.b;p.r1;p.r2;p.r3;p.a;p.d];
x0 = log(x0);
nip = length(x0);
for jj = 1:length(X(:))

    p.alpha = X(jj);
    p.beta  = Y(jj);
    %fprintf('alpha = %3.3f',p.alpha)
    L = @(x) neglogpost_log(x,p,grd,M2d);

    options = optimoptions(@fminunc,'Algorithm','trust-region',...
        'GradObj','on','Hessian','on','Display','off',...
        'MaxFunEvals',1000,'MaxIter',1000,'TolX',1e-7,...
        'DerivativeCheck','off','FinDiffType', ...
        'central','TolFun',1e-7,'PrecondBandWidth',Inf);

    [xhat,fval,exitflag] = fminunc(L,x0,options);

    [f,dfdx,d2fdx2] = neglogpost_log(xhat,p,grd,M2d);
    HH = d2fdx2;
    logZ(jj) = -f-0.5*log(det(HH))+6/2*log(p.alpha)+(6*idata)/2*log(p.beta);
end

%figure(2)
%contourf(X,Y,logZ);colorbar
%set(gca,'XTick',[1:1:40])
%set(gca,'XTickLabel',{'0','10','20','30','40','50','60','70','80','90','100'})
%set(gca,'YTick',[20:1:60])
%set(gca,'YTickLabel',{'0','10','20','30','40','50','60','70','80','90','100'})
%text('Interpreter','latex','String','$$\Lambda = 6.12$$','Position',[35 55])
%text('Interpreter','latex','String','$$\Gamma = 36.97$$','Position',[35 53])
%text(5,55,'MedFlux')
%xlabel('\Lambda-scaling factor for parameter')
%ylabel('\Gamma-scaling factor for data')

imax = find(logZ(:)==max(logZ(:)));
p.alpha = X(imax);
p.beta  = Y(imax);
%keyboard
[xhat,fval,exitflag] = fminunc(L,x0,options);
[f,dfdx,d2fdx2] = neglogpost_log(xhat,p,grd,M2d);
HH = d2fdx2;

I = speye(nip);
N_p = nip-alpha*(trace(HH.^(-1)*I));

error = sqrt(diag(inv(HH)));
R.logZ = logZ;
R.upbar = (exp(xhat+error)-exp(xhat));
R.lowbar = (exp(xhat)-exp(xhat-error));
R.xhat = exp(xhat);
R.alpha = p.alpha;
R.beta = p.beta;
fname = sprintf('xhat_log');
save(fname,'R');
