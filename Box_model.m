% this script is used to load data, set up constants and 
% parameters. It also sets up optimization conditions for 
% Matlab routine fminunc. At the end, it calculates errorbars
% for optimal parameter values and save them to file. 
addpath('~/Dropbox/MOCM/WEILEI/myfunc')
load DB2_May2005.mat
p.dp   = depth;     % sampling depths.
p.Chl  = Chla;      % large sized Chl a.
p.POC  = POC;       % large sized particulate organic matter.
p.Phyo = Phyeo;     % large sized phaeopigment.
p.chl  = chla;      % small sized Chl a.
p.poc  = poc;       % small sized particulated organic matter.
p.phyo = phyeo;     % small sized phaeopigment.

% calculate standard deviation for each data group. And
% they will be used to normalize fitting errors for each group.
p.Chl_std  = std(p.Chl);
p.chl_std  = std(p.chl);
p.POC_std  = std(p.POC);
p.poc_std  = std(p.poc);
p.Phyo_std = std(p.Phyo);
p.phyo_std = std(p.phyo);
idata      = length(p.Chl);
% second per year.
spy        = 365*24*60*60;

jj = length(p.dp);
M2d = ones(1,jj);

grd = buildgrd(p,M2d);
p.b   = 0.84;   % Martin Curve exponential
p.r1  = 1;      % small sized Chl a to phyeopigment rate constant [y^[-1]]
p.r2  = 1;      % small sized POC remineralization rate [y^[-1]];
p.r3  = 1;      % Phyeopigment reminearlization rate [y^[-1]];
p.a   = 3;      % aggregation rate      [y^[-1]];
p.d   = 150;    % disaggregation rate   [y^[-1]];

p.eta = 1.0;     % coefficient to convert conc. to production rate

load xhat_var_SV.mat

% setting up grid search for alpha and beta.
%alpha = linspace(1,5,40);
%beta  = linspace(15,35,20);
alpha = R.alpha;
beta  = R.beta;
[X,Y] = meshgrid(alpha,beta);

logZ = 0*X;
x0 = R.xhat;
%x0 = [b;r1;r2;r3;a;d];
x0 = log(x0);
nip = length(x0);
for jj = 1:length(X(:))

    p.alpha = X(jj);
    p.beta  = Y(jj);

    L = @(x) neglogpost(x,p,grd,M2d);

    options = optimoptions(@fminunc,'Algorithm','trust-region',...
        'GradObj','on','Hessian','on','Display','off',...
        'MaxFunEvals',1000,'MaxIter',1000,'TolX',1e-9,...
        'DerivativeCheck','off','FinDiffType', ...
        'central','TolFun',1e-9,'PrecondBandWidth',Inf);

    [xhat,fval,exitflag] = fminunc(L,x0,options);

    [f,dfdx,d2fdx2] = neglogpost(xhat,p,grd,M2d);
    HH = d2fdx2;
    logZ(jj) = -f-0.5*log(det(HH))+6/2*log(p.alpha)+(6*idata)/2*log(p.beta);
end

%%%%%%%% Comment this out to get a contour plot for alpha and beta %%%%%%%%%%%%
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
%%%%%%%%  Comment this out to get a contour plot for alpha and beta %%%%%%%%%%% 


% finding optimal alpha and beta, recalculate
% Hessian matrix based on them, and calculate error bars for
% each parameter values.
imax = find(logZ(:)==max(logZ(:)));
p.alpha = X(imax);
p.beta  = Y(imax);

[xhat,fval,exitflag] = fminunc(L,x0,options);
[f,dfdx,d2fdx2] = neglogpost(xhat,p,grd,M2d);
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
fname = sprintf('xhat_var_SV');
save(fname,'R');
