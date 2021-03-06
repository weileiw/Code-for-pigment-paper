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
p.dp   = depth;
p.Chl_std  = std(Chla);
p.POC_std  = std(POC);
p.Phyo_std = std(Phyeo);
p.chl_std  = std(chla);
p.poc_std  = std(poc);
p.phyo_std = std(phyeo);
% second per year
spy      = 365*24*60*60;

depth = p.dp;
jj = length(depth);
M2d = ones(1,jj);

grd = buildgrd(p,M2d);
p.b   = 0.84;   % Martin Curve exponential
p.r1  = 1;      % small sized Chl a to phyeopigment rate constant [y^[-1]]
p.r2  = 1;      % small sized POC remineralization rate [y^[-1]];
p.r3  = 1;      % Phyeopigment reminearlization rate [y^[-1]];
p.a   = 3;      % aggregation rate      [y^[-1]];
p.d   = 150;    % disaggregation rate   [y^[-1]];

p.eta = 1.0;      % coefficient to convert conc. to production rate [dimensionless];

load xhat_4p.mat

% setting up grid search for alpha and beta.
%alpha = linspace( 1,80,80); %range for testing eta;
%beta  = linspace(16,40,50); %range for testing eta;
%alpha = linspace(11,60,100); %range for eta = 1;
%beta  = linspace(16,40,50); %range for eta = 1;
alpha = R.alpha;
beta  = R.beta;

[X,Y] = meshgrid(alpha,beta);

logZ = 0*X;
x0 = R.xhat;
%x0  = [p.b;p.r1;p.a;p.b]; 
x0  = log(x0);
nip = length(x0);

idata = length(depth)-1;
for jj = 1:length(X(:))

    p.alpha = X(jj);
    p.beta  = Y(jj);
    L = @(x) neglogpost_4p(x,p,grd,M2d);

    options = optimoptions(@fminunc,'Algorithm','trust-region',...
        'GradObj','on','Hessian','on','Display','off',...
        'MaxFunEvals',1000,'MaxIter',1000,'TolX',1e-9,...
        'DerivativeCheck','off','FinDiffType', ...
        'central','TolFun',1e-9,'PrecondBandWidth',Inf);

    [xhat,fval,exitflag] = fminunc(L,x0,options);

    [f,dfdx,d2fdx2] = neglogpost_4p(xhat,p,grd,M2d);
    HH = d2fdx2;
    logZ(jj) = -f-0.5*log(det(HH))+nip/2*log(p.alpha)+(6*idata)/2*log(p.beta);
end

%%%%%%%% Comment this out to get a contour plot for alpha and beta %%%%%%%%%%%%
%figure(2)
contourf(X,Y,logZ);colorbar
%set(gca,'XTick',[0:1:10])
%set(gca,'XTickLabel',{'0','1','2','3','4','5','6','7','8','9','10'})
%set(gca,'YTick',[20:5:60])
%set(gca,'YTickLabel',{'20','25','30','35','40','45','50','55','60'})
%text('Interpreter','latex','String','$$\Lambda = 21.89$$','Position',[50 35],'fontsize',16)
%text('Interpreter','latex','String','$$\Gamma = 27.76$$','Position',[50 37],'fontsize',16)
%text(5,55,'MedFlux')
xlabel('\Lambda-scaling factor for parameter')
ylabel('\Gamma-scaling factor for data')
set(gca,'fontsize',16)
%%%%%%%%  Comment this out to get a contour plot for alpha and beta %%%%%%%%%%% 

% finding optimal alpha and beta, recalculate
% Hessian matrix based on them, and calculate error bars for
% each parameter values.
imax = find(logZ(:)==max(logZ(:)));
p.alpha = X(imax);
p.beta  = Y(imax);
% recalculate xhat based on optimized alpha and beta.
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
fname = sprintf('xhat_4p');
%save(fname,'R');
