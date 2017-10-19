clc
clear all
%
%
Dab = 211; Dac = 1605; Dbc = 1394;


iteration = 450000;
changefact = 0.01;
tempint = 1;
printint = 1000;
relimitint = 10000;

% ===========================Reading parameters in, the limit of each
% paremeters were initially set in a arbitray way, and were changed
% accordingly based on the first run results.

[pdata, ptext] = xlsread('parameter pigment model.xlsx');
pname = ptext(:,1);
pvalue = pdata(:,1);
pmin = pdata(:,3);
pmax = pdata(:,4);


for i = 1:length(pname)
    
    if strcmp (pname(i), 'k1a')
        
        k1a = pvalue(i);
        mink1a = pmin(i);
        maxk1a = pmax(i);
        
    elseif strcmp (pname(i),'k2a')
        
        k2a = pvalue(i);
        mink2a = pmin(i);
        maxk2a = pmax(i);
        
    elseif strcmp (pname(i),'t1a')
        t1a = pvalue(i);
        mint1a = pmin(i);
        maxt1a = pmax(i);
        
    elseif strcmp (pname(i),'t2a')
        
        t2a = pvalue(i);
        mint2a = pmin(i);
        maxt2a = pmax(i);
        
    elseif strcmp (pname(i),'k1b')
        
        k1b = pvalue(i);
        mink1b = pmin(i);
        maxk1b = pmax(i);
        
    elseif strcmp (pname(i),'k2b')
        
        k2b = pvalue(i);
        mink2b = pmin(i);
        maxk2b = pmax(i);
        
    elseif strcmp (pname(i),'t1b')
        t1b = pvalue(i);
        mint1b = pmin(i);
        maxt1b = pmax(i);
        
    elseif strcmp (pname(i),'t2b')
        
        t2b = pvalue(i);
        mint2b = pmin(i);
        maxt2b = pmax(i);
        
        
    end
    
end



ptemp=ones(1,12)*tempint;

% =====================Reading data in.

filename='Pigment data in matlab';
datasheet='49';

[ddata, dtext] = xlsread(filename,datasheet);
dname = dtext((2:length(dtext)),1);
dvalue = ddata(:,1);

for i = 1:length(dname)
    
    
    if strcmp(dname(i),'Chla')
        Chla = dvalue(i);
        
    elseif strcmp(dname(i),'ChlA')
        ChlA = dvalue(i);
        
    elseif strcmp(dname(i),'Pra')
        Pra = dvalue(i);
        
    elseif strcmp(dname(i),'PrA')
        PrA = dvalue(i);
        
    elseif strcmp(dname(i),'Chlb')
        Chlb = dvalue(i);
        
    elseif strcmp(dname(i),'ChlB')
        ChlB = dvalue(i);
        
    elseif strcmp(dname(i),'Prb')
        Prb = dvalue(i);
        
    elseif strcmp(dname(i),'PrB')
        PrB = dvalue(i);
        
    elseif strcmp(dname(i),'Chlc')
        Chlc = dvalue(i);
        
    elseif strcmp(dname(i),'ChlC')
        ChlC = dvalue(i);
        
    elseif strcmp(dname(i),'Prc')
        Prc = dvalue(i);
        
    elseif strcmp(dname(i),'PrC')
        PrC = dvalue(i);
        
    elseif strcmp(dname(i),'Pa')
        Pa = dvalue(i);
        
    elseif strcmp(dname(i),'PA')
        PA = dvalue(i);
        
    elseif strcmp(dname(i),'Pb')
        Pb = dvalue(i);
        
    elseif strcmp(dname(i),'PB')
        PB = dvalue(i);
        
    elseif strcmp(dname(i),'Pc')
        Pc = dvalue(i);
        
    elseif strcmp(dname(i),'PC')
        PC = dvalue(i);
        
    elseif strcmp(dname(i),'Fa')
        Fa = dvalue(i);
        
    elseif strcmp(dname(i),'FA')      
        FA = dvalue(i);
        
    elseif strcmp(dname(i),'Fb')
        Fb = dvalue(i);
        
    elseif strcmp(dname(i),'FB')
        FB = dvalue(i);
        
    elseif strcmp(dname(i),'Fc')
        Fc = dvalue(i);
        
    elseif strcmp(dname(i),'FC')
        FC = dvalue(i);
        
    elseif strcmp(dname(i), 'Fchla')
        Fchla=dvalue(i);
        
    elseif strcmp(dname(i), 'FchlA')
        FchlA=dvalue(i);
        
    elseif strcmp(dname(i), 'Fpra')
        Fpra=dvalue(i);
        
    elseif strcmp(dname(i), 'FprA')
        FprA=dvalue(i);
        
    elseif strcmp(dname(i), 'Fchlb')
        Fchlb=dvalue(i);
        
    elseif strcmp(dname(i), 'FchlB')
        FchlB=dvalue(i);
        
    elseif strcmp(dname(i), 'Fprb')
        Fprb=dvalue(i);
        
    elseif strcmp(dname(i), 'FprB')
        FprB=dvalue(i);
        
    elseif strcmp(dname(i), 'Fchlc')
        Fchlc=dvalue(i);
        
    elseif strcmp(dname(i), 'FchlC')
        FchlC=dvalue(i);
        
    elseif strcmp(dname(i), 'Fprc')
        Fprc=dvalue(i);
        
    elseif strcmp(dname(i), 'FprC')
        FprC=dvalue(i);
        
    end
    
end

p = [k1a,k2a,t1a,t2a,k1b,k2b,t1b,t2b];
oldp = [k1a,k2a,t1a,t2a,k1b,k2b,t1b,t2b];
maxp = [k1a,k2a,t1a,t2a,k1b,k2b,t1b,t2b];

prmax = [maxk1a,maxk2a,maxt1a,maxt2a,maxk1b,maxk2b,maxt1b,maxt2b];
prmin = [mink1a,mink2a,mint1a,mint2a,mink1b,mink2b,mint1b,mint2b];

paramptr = 1;
numparams = 8;
za = 313; zb = 524; zc = 1918;

ChlA = 0.5*(zb*(ChlA+ChlB)+za*(ChlB-3*ChlA))/(zb-za);
ChlB = 0.5*(zc*(ChlB+ChlC)+zb*(ChlC-3*ChlB))/(zc-zb);
Chla = 0.5*(zb*(Chla+Chlb)+za*(Chlb-3*Chla))/(zb-za);
Chlb = 0.5*(zc*(Chlb+Chlc)+zb*(Chlc-3*Chlb))/(zc-zb);

PR1 = 0.5*(zb*(PrA+PrB)+za*(PrB-3*PrA))/(zb-za);
PR2 = 0.5*(zc*(PrB+PrC)+zb*(PrC-3*PrB))/(zc-zb);
pr1 = 0.5*(zb*(Pra+Prb)+za*(Prb-3*Pra))/(zb-za);
pr2 = 0.5*(zc*(Prb+Prc)+zb*(Prc-3*Prb))/(zc-zb);

P1 = 0.5*(zb*(PA+PB)+za*(PB-3*PA))/(zb-za);
P2 = 0.5*(zc*(PB+PC)+zb*(PC-3*PB))/(zc-zb);
p1 = 0.5*(zb*(Pa+Pb)+za*(Pb-3*Pa))/(zb-za);
p2 = 0.5*(zc*(Pb+Pc)+zb*(Pc-3*Pb))/(zc-zb);


Fla = (Fchlb-Fchla)/Dab;
Flb = (Fchlc-Fchlb)/Dbc;

Fra = (Fprb-Fpra)/Dab;
Frb = (Fprc-Fprb)/Dbc;

Fpa = (Fb-Fa)/Dab;
Fpb = (Fc-Fb)/Dbc;

FlA = (FchlB-FchlA)/Dab;
FlB = (FchlC-FchlB)/Dbc;

FrA = (FprB-FprA)/Dab;
FrB = (FprC-FprB)/Dbc;

FpA = (FB-FA)/Dab;
FpB = (FC-FB)/Dbc;


for iter = 1:iteration
    
    % k1 aggregation rate constant;
    % k2 disaggregation rate constant;
    % t1 chl a degradation rate constant;
    % t2 OM degradation rate constant;
    p = [k1a,k2a,t1a,t2a,k1b,k2b,t1b,t2b];
    

    eFla = ChlA*k2a-Chla*(k1a+t1a);
    eFlb = ChlB*k2b-Chlb*(k1b+t1b);
    
    eFra = Chla*t1a+PR1*k2a-pr1*(t2a+k1a);
    eFrb = Chlb*t1b+PR2*k2b-pr2*(t2b+k1b);
    
    eFpa = P1*k2a-p1*(k1a+t2a);
    eFpb = P2*k2b-p2*(k1b+t2b);
    
    eFlA = Chla*k1a-ChlA*k2a;
    eFlB = Chlb*k1b-ChlB*k2b;
    
    eFrA = pr1*k1a-PR1*k2a;
    eFrB = pr2*k1b-PR2*k2b;
    
    eFpA = p1*k1a-P1*k2a;
    eFpB = p2*k1b-P2*k2b;
    
%     measurementl = [Fla Flb Flc FlA FlB FlC];
%     measurementr = [Fra Frb Frc FrA FrB FrC];
%     measurementp = [Fpa Fpb Fpc FpA FpB FpC];
%     
%     estimationl = [eFla eFlb eFlc eFlA eFlB eFlC];
%     estimationr = [eFra eFrb eFrc eFrA eFrB eFrC];
%     estimationp = [eFpa eFpb eFpc eFpA eFpB eFpC];
    
    measurements = [Fla Flb FlA FlB Fra Frb FrA FrB Fpa Fpb FpA FpB];
    estimations = [eFla eFlb eFlA eFlB eFra eFrb eFrA eFrB eFpa eFpb eFpA eFpB];
    
    L = length(measurements);
    loglike = -(L/2)*log(sum((measurements-estimations).^2)/L);
    
%         subloglike(1) = -(L1/2)*log(sum((measurementl-estimationl).^2)/L1);
%         subloglike(2) = -(L2/2)*log(sum((measurementr-estimationr).^2)/L2);
%         subloglike(3) = -(L3/2)*log(sum((measurementp-estimationp).^2)/L3);
%         
%         loglike = sum (subloglike);
    

    if loglike == Inf
        
        break;
        
    end
    
    
    
    if iter == 1
        
        oldloglike = loglike;
        maxloglike = loglike;
        accept='n';
        
        
    else
        
        
        logdiff = loglike-oldloglike;
        accept='n';
        
        
        if logdiff>0
            accept = 'y';
            
            % elseif (logdiff>-50 && exp(logdiff/ptemp(paramptr))>rand)
        elseif exp(logdiff/ptemp(paramptr))>rand
            accept = 'y';
            
        end
    end
    
    
    if accept == 'y' %&& flag=='y';
        
        oldloglike = loglike;
        oldp(paramptr) = p(paramptr);
        
        %Armstrong
        % ptemp(paramptr)=ptemp(paramptr)/1.01;
        %
        ptemp(paramptr) = ptemp(paramptr)/1.01;
        
        
        if loglike > maxloglike
            
            maxloglike = loglike;
            
            for k = 1:numparams
                
                maxp(k) = p(k);
                
            end
        end
        
    else
        
        p(paramptr) = oldp(paramptr);
        ptemp(paramptr) = ptemp(paramptr)*1.01;
        
    end
    
    paramptr = randi(numparams);
    
    
    
    while 0<1
        
        % p(paramptr)=roundn(oldp(paramptr)*exp(changefact*log(prmax(paramptr)/prmin(paramptr))*(rand-0.5)),-9);
        
        p(paramptr) = oldp(paramptr)*exp(changefact*log(prmax(paramptr)/prmin(paramptr))*(rand-0.5));
        
        i = i+1;
        % p(paramptr)=oldp(paramptr)*exp(changefact*ptemp(paramptr)*(rand-0.5));
        
        if p(paramptr) < prmax(paramptr) && p(paramptr) > prmin(paramptr)
            
            break;
        end
        
        
    end
    
    k1a=p(1); k2a=p(2); t1a=p(3); t2a=p(4); k1b=p(5); k2b=p(6); t1b=p(7); t2b=p(8);
     
    
    if (rem(iter, printint)==0) && (iter~=0)
        
        for i=1:numparams
            
            disp([pname(i),p(i),maxp(i),ptemp(i)])
            
        end
        
        X2=sprintf('loglike is %d',loglike);
        X3=sprintf('maxloglike is %d', maxloglike);
        disp(X2)
        disp(X3)
        
    end
    
    
    if (rem(iter, relimitint)==0)&&(iter~=0)&&(iter~=iteration)
        
        for i=1:numparams
            
            p(i)=maxp(i);
            oldp(i)=maxp(i);
            prmin(i)=sqrt(prmin(i)*maxp(i));
            prmax(i)=sqrt(prmax(i)*maxp(i));
            
        end
        
        %reassign the parameters, and bring this parameters into the
        %calculations.
        
        k1a=p(1); k2a=p(2); t1a=p(3); t2a=p(4); k1b=p(5); k2b=p(6); t1b=p(7); t2b=p(8);
 
        
        Y1=sprintf('number of iteration is: %d', iter);
        Y2=sprintf('current likelyhood is: %d',loglike);
        Y3=sprintf('maxlikelyhood is: %d', maxloglike);
        disp(Y1)
        disp(Y2)
        disp(Y3)
        
        
    end
end

% Z1=sprintf('loglike is %d and maxloglike is %d:', loglike, maxloglike);
% disp(Z1)
% display(p)

%
%


name={'maxloglike','loglike','k1a','k2a','t1a','t2a','k1b','k2b','t1b','t2b'};

R=[maxloglike, loglike, maxp];
N=cell(size(name'));
V=cell(size(maxp));

for i=1:numel(name)
    
    N{i}=char(name(i));
    V{i}=R(i);
    
end


if exist ('opt.mat') == 0
    out = [];
else 
    load opt.mat out
    out = out;
end

opt = [k1a,k2a,t1a,t2a,...
       k1b,k2b,t1b,t2b];
out = cat(1,out,opt);

fname = sprintf('opt');
save(fname,'out');
