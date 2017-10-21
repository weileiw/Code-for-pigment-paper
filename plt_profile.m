load DB2_May2005.mat 

figure()
plot(phyeo, depth,'b*')
hold on 
plot(chla, depth,'r+')
hold off 
legend('small-sized phyeo','small-sized Chla') 
xlabel('Phyeo/Chla concentration (\muM)','fontsize',16)
ylabel('Dpeth (m)','fontsize', 16)
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
set(gca,'fontsize',16)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
print -depsc chl_phyeo_profile.eps

figure()
plot(Phyeo, depth,'b*')
hold on 
plot(Chla, depth,'r+')
hold off 
legend('large-sized Phyeo','large-sized Chla') 
xlabel('Phyeo/Chla concentration (\muM)','fontsize',16)
ylabel('Dpeth (m)','fontsize', 16)
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
set(gca,'fontsize',16)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
print -depsc Phyeo_Chla_profile.eps

figure()
plot(POC, depth,'b*')
hold on 
plot(poc, depth,'r+')
hold off 
legend('large-sized POC','small-sized POC') 
xlabel('POC concentration (\muM)','fontsize',16)
ylabel('Dpeth (m)','fontsize', 16)
set(gca,'XAxisLocation','top','YAxisLocation','left','ydir','reverse');
set(gca,'fontsize',16)
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
print -depsc POC_poc_profile.eps

