clear all,close all ,clc
b=[0.008
0.0111
0.0089
0.01
0.006
0.007872362
];
r=[0.195056809
0.1841
0.1901
0.1708
0.2214
0.179948292];
uv1=[0.1 0.25 0.5 0.75 0];
T=0.25;
t=0:0.25:12;
N = size(t,2);
S=zeros(length(b),length(t),length(uv1));
I=zeros(length(b),length(N),length(uv1));
z=zeros(length(b),length(N),length(uv1));
for m=1:1:length(uv1)
    S(:,1,m)=[70;66;70;76;67;77]; %[x_WP; x_MIA; x_FTM; x_TP; x_MCO; x_JAX];
    I(:,1,m)=[30;34;30;24;33;23];
    z(:,1,m)=[30*uv1(m);34*uv1(m);30*uv1(m);24*uv1(m);33*uv1(m);23*uv1(m)];
    z1(:,1,m)=[0.30*uv1(m);0.34*uv1(m);0.30*uv1(m);0.24*uv1(m);0.33*uv1(m);0.23*uv1(m)];
end


for j=1:1:length(uv1)
for i=1:1:length(b)
for k = 2:N
    if S(i,k-1,j) >= (r(i)/b(i))
        uv(i,k-1,j)=uv1(j);
    else
        uv(i,k-1,j)=0;
    end

    S(i,k,j)= S(i,k-1,j)+((-b(i)*S(i,k-1,j)*I(i,k-1,j)-uv(i,k-1,j)*S(i,k-1,j))*T);
    I(i,k,j)= I(i,k-1,j)+(b(i)*S(i,k-1,j)*I(i,k-1,j)-r(i)*I(i,k-1,j))*T;
    z(i,k,j)=(z(i,k-1,j)+((uv(i,k-1,j)*S(i,k-1,j))*T));
    z1(i,k,j)=(z1(i,k-1,j)+((uv(i,k-1,j)*15e3*831*0.01*S(i,k-1,j))*T));
    
end
end
end
ts=[];
Imax=[];
for j=1:1:length(uv1)
for i=1:1:length(b)
    Imax(i,j)=max(I(i,:,j)); 

end
end

for j=1:1:length(uv1)
for i=1:1:length(b)
for k = 1:N-1
    if uv(i,k,j) == 0
       ts(i,j)=t(k)-T;
    break;
    end
    
end
end
end
time=0:0.001:12;
sab=size(0:0.001:12,2);
uvnew=zeros(sab,length(uv1),length(b));

for j=1:1:length(uv1)
for i=1:1:length(b)
a=size(0:0.001:ts(i,j),2);
uvnew(1:a,i,j)=uv1(j);
end
end



city={'West Palm' 'Miami/FT Lauderdale' 'Fort Myers Naples' 'Tampa St Pete' 'Orlando' 'Jacksonville'};

% for c=1:1:length(b)
% for a=1:1:length(uv1)
%     figure(c)
%     subplot(3,1,1)
%     plot(t,S(c,:,a))
%     xlabel('Days')
%     ylabel('S(t) % of Gas Stations')
%     title(sprintf('Optimal Refuelling Policy for %s', city{c}))
%     hold on
%     subplot(3,1,3)
%     plot(t,I(c,:,a))
%     xlabel('Days')
%     ylabel('I(t) % of Gas Stations')
%     hold on
%     subplot (3,1,2)
%     plot(time(1:end),uvnew(:,c,a))
%     legend ('uv=0.1','uv=0.25','uv=0.5','uv=0.75','uv=0')
%     xlabel('Days')
%     hold on
% %     subplot(4,1,4)
% %     plot(t,z(c,:,a))
% %     ylabel('Cost')
% %     xlabel('Days')
% %     hold on
% end
% end
% clc
% 
% figure
p=2;
% plot(t,S(p,:,5),'-k+',t,S(p,:,1),'-ro',t,S(p,:,2),'-bx',t,S(p,:,3),'-g^',t,S(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 5])
% ylabel('S(t) [%]')
% xlabel('Days')
% sheet='MiamiFt Lauderdale Area';
% fn = sprintf('Irma_OC_S%s',sheet); 
% saveas( gcf, fn, 'png' );

% 
% 
% figure
% plot(t,I(2,:,5),'-k+',t,I(2,:,1),'-ro',t,I(2,:,2),'-bx',t,I(2,:,3),'-g^',t,I(2,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 10])
% ylabel('I(t) [%]')
% xlabel('Days')
% sheet='MiamiFt Lauderdale Area';
% fn = sprintf('Irma_OC_I_%s',sheet); 
% saveas( gcf, fn, 'png' );
% 
% 
% 
% figure
% plot(time,uvnew(:,p,5),'-k',time,uvnew(:,p,1),'-r',time,uvnew(:,p,2),'-b',time,uvnew(:,p,3),'-g',time,uvnew(:,p,4),'-m')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% title('u_{v,max} for Miami/Fort Lauderdale')
% xlabel('Time(s)')
% ylabel('u_v (t)')
% xlim([0 10])
q=5;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MIA_OC_Data.xls', e, 'uv=0', 'A2')


q=1;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MIA_OC_Data.xls', e, 'uv=0.1', 'A2')
xlswrite('MIA_OC_Data.xls', d(q), 'uv=0.1', 'D2')

q=2;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MIA_OC_Data.xls', e, 'uv=0.25', 'A2')
xlswrite('MIA_OC_Data.xls', d(q), 'uv=0.25', 'D2')

q=3;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MIA_OC_Data.xls', e, 'uv=0.5', 'A2')
xlswrite('MIA_OC_Data.xls', d(q), 'uv=0.5', 'D2')

q=4;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MIA_OC_Data.xls', e, 'uv=0.75', 'A2')
xlswrite('MIA_OC_Data.xls', d(q), 'uv=0.75', 'D2')


% 
% city={'West Palm' 'Miami/FT Lauderdale' 'Fort Myers Naples' 'Tampa St Pete' 'Orlando' 'Jacksonville'};
% %West Palm
 p=1;
 q=5;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('WP_OC_Data.xls', e, 'uv=0', 'A2')


q=1;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('WP_OC_Data.xls', e, 'uv=0.1', 'A2')
xlswrite('WP_OC_Data.xls', d(q), 'uv=0.1', 'D2')

q=2;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('WP_OC_Data.xls', e, 'uv=0.25', 'A2')
xlswrite('WP_OC_Data.xls', d(q), 'uv=0.25', 'D2')

q=3;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('WP_OC_Data.xls', e, 'uv=0.5', 'A2')
xlswrite('WP_OC_Data.xls', d(q), 'uv=0.5', 'D2')

q=4;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('WP_OC_Data.xls', e, 'uv=0.75', 'A2')
xlswrite('WP_OC_Data.xls', d(q), 'uv=0.75', 'D2')
% figure
% plot(t,S(p,:,5),'-k+',t,S(p,:,1),'-ro',t,S(p,:,2),'-bx',t,S(p,:,3),'-g^',t,S(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 5])
% xlabel('Time(s)')
% ylabel('S(t) [%]')
% sheet='West Palm';
% fn = sprintf('Irma_OC_S_%s',sheet); 
% saveas( gcf, fn, 'png' );
% 
% 
% figure
% plot(t,I(p,:,5),'-k+',t,I(p,:,1),'-ro',t,I(p,:,2),'-bx',t,I(p,:,3),'-g^',t,I(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 10])
% ylabel('I(t) [%]')
% xlabel('Days')
% sheet='West Palm';
% fn = sprintf('Irma_OC_I_%s',sheet); 
% saveas( gcf, fn, 'png' );
% 
% 
% figure
% plot(time,uvnew(:,p,5),'-k',time,uvnew(:,p,1),'-r',time,uvnew(:,p,2),'-b',time,uvnew(:,p,3),'-g',time,uvnew(:,p,4),'-m')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% title('u_{v,max} for West Palm')
% xlabel('Days')
% ylabel('u_v (t)')
% xlim([0 10])
% 
% %Fort Myers Naples
 p=3;
 q=5;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('FTM_OC_Data.xls', e, 'uv=0', 'A2')


q=1;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('FTM_OC_Data.xls', e, 'uv=0.1', 'A2')
xlswrite('FTM_OC_Data.xls', d(q), 'uv=0.1', 'D2')

q=2;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('FTM_OC_Data.xls', e, 'uv=0.25', 'A2')
xlswrite('FTM_OC_Data.xls', d(q), 'uv=0.25', 'D2')

q=3;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('FTM_OC_Data.xls', e, 'uv=0.5', 'A2')
xlswrite('FTM_OC_Data.xls', d(q), 'uv=0.5', 'D2')

q=4;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('FTM_OC_Data.xls', e, 'uv=0.75', 'A2')
xlswrite('FTM_OC_Data.xls', d(q), 'uv=0.75', 'D2')
% sheet = 'FortMyersNaples';
% figure
% plot(t,S(p,:,5),'-k+',t,S(p,:,1),'-ro',t,S(p,:,2),'-bx',t,S(p,:,3),'-g^',t,S(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 5])
% xlabel('Days')
% ylabel('S(t) [%]')
% fn = sprintf('Irma_OC_S_%s',sheet); 
% saveas( gcf, fn, 'png' );
% 
% 
% figure
% plot(t,I(p,:,5),'-k+',t,I(p,:,1),'-ro',t,I(p,:,2),'-bx',t,I(p,:,3),'-g^',t,I(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 10])
% ylabel('I(t) [%]')
% xlabel('Days')
% fn = sprintf('Irma_OC_I_%s',sheet); 
% saveas( gcf, fn, 'png' );
% 
% 
% figure
% plot(time,uvnew(:,p,5),'-k',time,uvnew(:,p,1),'-r',time,uvnew(:,p,2),'-b',time,uvnew(:,p,3),'-g',time,uvnew(:,p,4),'-m')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlabel('Days')
% ylabel('u_v (t)')
% xlim([0 10])
% fn = sprintf('Irma_OC_uv_%s',sheet); 
% saveas( gcf, fn, 'png' );
% 
% %Tampa/St Pete
 p=4;
  q=5;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('TPA_OC_Data.xls', e, 'uv=0', 'A2')


q=1;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('TPA_OC_Data.xls', e, 'uv=0.1', 'A2')
xlswrite('TPA_OC_Data.xls', d(q), 'uv=0.1', 'D2')

q=2;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('TPA_OC_Data.xls', e, 'uv=0.25', 'A2')
xlswrite('TPA_OC_Data.xls', d(q), 'uv=0.25', 'D2')

q=3;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('TPA_OC_Data.xls', e, 'uv=0.5', 'A2')
xlswrite('TPA_OC_Data.xls', d(q), 'uv=0.5', 'D2')

q=4;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('TPA_OC_Data.xls', e, 'uv=0.75', 'A2')
xlswrite('TPA_OC_Data.xls', d(q), 'uv=0.75', 'D2')
% figure
% plot(t,S(p,:,5),'-k+',t,S(p,:,1),'-ro',t,S(p,:,2),'-bx',t,S(p,:,3),'-g^',t,S(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 5])
% xlabel('Days')
% ylabel('S(t) [%]')
% sheet='TampaStPete';
% fn = sprintf('Irma_OC_S_%s',sheet); 
% saveas( gcf, fn, 'png' )
% 
% 
% 
% figure
% plot(t,I(p,:,5),'-k+',t,I(p,:,1),'-ro',t,I(p,:,2),'-bx',t,I(p,:,3),'-g^',t,I(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 10])
% ylabel('I(t) [%]')
% xlabel('Days')
% fn = sprintf('Irma_OC_I_%s',sheet); 
% saveas( gcf, fn, 'png' )
% 
% 
% figure
% plot(time,uvnew(:,p,5),'-k',time,uvnew(:,p,1),'-r',time,uvnew(:,p,2),'-b',time,uvnew(:,p,3),'-g',time,uvnew(:,p,4),'-m')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% title('u_{v,max} for Tampa/St Pete')
% xlabel('Time(s)')
% ylabel('u_v (t)')
% xlim([0 10])
% 
% %Orlando
 p=5;
  q=5;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MCO_OC_Data.xls', e, 'uv=0', 'A2')


q=1;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MCO_OC_Data.xls', e, 'uv=0.1', 'A2')
xlswrite('MCO_OC_Data.xls', d(q), 'uv=0.1', 'D2')

q=2;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MCO_OC_Data.xls', e, 'uv=0.25', 'A2')
xlswrite('MCO_OC_Data.xls', d(q), 'uv=0.25', 'D2')

q=3;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MCO_OC_Data.xls', e, 'uv=0.5', 'A2')
xlswrite('MCO_OC_Data.xls', d(q), 'uv=0.5', 'D2')

q=4;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('MCO_OC_Data.xls', e, 'uv=0.75', 'A2')
xlswrite('MCO_OC_Data.xls', d(q), 'uv=0.75', 'D2')
% figure
% plot(t,S(p,:,5),'-k+',t,S(p,:,1),'-ro',t,S(p,:,2),'-bx',t,S(p,:,3),'-g^',t,S(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 5])
% xlabel('Days')
% ylabel('S(t) [%]')
% sheet='Orlando';
% fn = sprintf('Irma_OC_S_%s',sheet); 
% saveas( gcf, fn, 'png') 
% 
% 
% figure
% plot(t,I(p,:,5),'-k+',t,I(p,:,1),'-ro',t,I(p,:,2),'-bx',t,I(p,:,3),'-g^',t,I(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 10])
% ylabel('I(t) [%]')
% xlabel('Days')
% fn = sprintf('Irma_OC_I_%s',sheet); 
% saveas( gcf, fn, 'png' )
% 
% 
% figure
% plot(time,uvnew(:,p,5),'-k',time,uvnew(:,p,1),'-r',time,uvnew(:,p,2),'-b',time,uvnew(:,p,3),'-g',time,uvnew(:,p,4),'-m')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% title('u_{v,max} for Orlando')
% xlabel('Days')
% ylabel('u_v (t)')
% xlim([0 10])
% 
% figure
% plot(t,z1(p,:,5),'-k+',t,z1(p,:,1),'-ro',t,z1(p,:,2),'-bx',t,z1(p,:,3),'-g^',t,z1(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% title('u_{v,max} for Orlando')
% xlabel('Days')
% ylabel('Gallons of Fuel')
% 
% 
% 
% %Jacksonville
 p=6;
 q=5;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('JAX_OC_Data.xls', e, 'uv=0', 'A2')


q=1;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('JAX_OC_Data.xls', e, 'uv=0.1', 'A2')
xlswrite('JAX_OC_Data.xls', d(q), 'uv=0.1', 'D2')

q=2;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('JAX_OC_Data.xls', e, 'uv=0.25', 'A2')
xlswrite('JAX_OC_Data.xls', d(q), 'uv=0.25', 'D2')

q=3;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('JAX_OC_Data.xls', e, 'uv=0.5', 'A2')
xlswrite('JAX_OC_Data.xls', d(q), 'uv=0.5', 'D2')

q=4;
e=[t',I(p,:,q)',S(p,:,q)'];
d=[ts(2,1) ts(2,2) ts(2,3) ts(2,4) ts(2,5)];
xlswrite('JAX_OC_Data.xls', e, 'uv=0.75', 'A2')
xlswrite('JAX_OC_Data.xls', d(q), 'uv=0.75', 'D2')
% figure
% plot(t,S(p,:,5),'-k+',t,S(p,:,1),'-ro',t,S(p,:,2),'-bx',t,S(p,:,3),'-g^',t,S(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 5])
% xlabel('Days')
% ylabel('S(t) [%]')
% sheet='Jax';
% fn = sprintf('Irma_OC_S_%s',sheet); 
% saveas( gcf, fn, 'png') 
% 
% figure
% plot(t,I(p,:,5),'-k+',t,I(p,:,1),'-ro',t,I(p,:,2),'-bx',t,I(p,:,3),'-g^',t,I(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 10])
% ylabel('I(t) [%]')
% xlabel('Days')
% fn = sprintf('Irma_OC_I_%s',sheet); 
% saveas( gcf, fn, 'png' )
% 
% 
% figure
% plot(time,uvnew(:,p,5),'-k',time,uvnew(:,p,1),'-r',time,uvnew(:,p,2),'-b',time,uvnew(:,p,3),'-g',time,uvnew(:,p,4),'-m')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% title('u_{v,max} for Jacksonville')
% xlabel('Days')
% ylabel('u_v (t)')
% xlim([0 10])

clc
b_f=[0.012 0.0143];
r_f=[0.0953 0.1543];
uv1_f=0:0.05:1;
T_f=1/24;
t_f=0:T_f:18.4583;
N_f = size(t_f,2);
S_f=zeros(length(b_f),length(t_f),length(uv1_f));
I_f=zeros(length(b_f),length(N_f),length(uv1_f));
z_f=zeros(length(b_f),length(N_f),length(uv1_f));
for m=1:1:length(uv1_f)
    S_f(:,1,m)=[92.08;96.10]; %[x_Wilm; x_GNW;];
    I_f(:,1,m)=[100-92.08;100-96.10];
%     z(:,1,m)=[30*uv1(m);34*uv1(m);30*uv1(m);24*uv1(m);33*uv1(m);23*uv1(m)];
%     z1(:,1,m)=[0.30*uv1(m);0.34*uv1(m);0.30*uv1(m);0.24*uv1(m);0.33*uv1(m);0.23*uv1(m)];
end


for j=1:1:length(uv1_f)
for i=1:1:length(b_f)
for k = 2:N_f
    if S_f(i,k-1,j) >= (r_f(i)/b_f(i))
        uv_f(i,k-1,j)=uv1_f(j);
    else
        uv_f(i,k-1,j)=0;
    end

    S_f(i,k,j)= S_f(i,k-1,j)+((-b_f(i)*S_f(i,k-1,j)*I_f(i,k-1,j)-uv_f(i,k-1,j)*S_f(i,k-1,j))*T_f);
    I_f(i,k,j)= I_f(i,k-1,j)+(b_f(i)*S_f(i,k-1,j)*I_f(i,k-1,j)-r_f(i)*I_f(i,k-1,j))*T_f;
%     z(i,k,j)=(z(i,k-1,j)+((uv(i,k-1,j)*S(i,k-1,j))*T));
%     z1(i,k,j)=(z1(i,k-1,j)+((uv(i,k-1,j)*15e3*831*0.01*S(i,k-1,j))*T));
    
end
end
end
ts_f=[];
Imax_f=[];
for j=1:1:length(uv1_f)
for i=1:1:length(b_f)
    Imax_f(i,j)=max(I_f(i,:,j)); 

end
end

for j=1:1:length(uv1_f)
for i=1:1:length(b_f)
for k = 1:N_f-1
    if uv_f(i,k,j) == 0
       ts_f(i,j)=t_f(k)-T_f;
    break;
    end
    
end
end
end
time_f=0:0.001:t_f(end);
sab_f=size(0:0.001:t_f(end),2);
uvnew_f=zeros(sab_f,length(uv1_f),length(b_f));

for j=1:1:length(uv1_f)
for i=1:1:length(b_f)
a=size(0:0.001:ts_f(i,j),2);
uvnew_f(1:a,i,j)=uv1_f(j);
end
end

% %Wilmington,NC
% p=1;
% sheet='Wilmington (Florence)';
% figure
% plot(t_f,S_f(p,:,5),'-k+',t_f,S_f(p,:,1),'-ro',t_f,S_f(p,:,2),'-bx',t_f,S_f(p,:,3),'-g^',t_f,S_f(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 5])
% xlabel('Days')
% ylabel('S(t) [%]')
% fn = sprintf('Irma_OC_S_%s',sheet); 
% saveas( gcf, fn, 'png' );
% 
% figure
% plot(t_f,I_f(p,:,5),'-k+',t_f,I_f(p,:,1),'-ro',t_f,I_f(p,:,2),'-bx',t_f,I_f(p,:,3),'-g^',t_f,I_f(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 10])
% ylabel('I(t) [%]')
% xlabel('Days')
% fn = sprintf('Irma_OC_I_%s',sheet); 
% saveas( gcf, fn, 'png' )
% 
% figure
% plot(time_f,uvnew_f(:,p,5),'-k',time_f,uvnew_f(:,p,1),'-r',time_f,uvnew_f(:,p,2),'-b',time_f,uvnew_f(:,p,3),'-g',time_f,uvnew_f(:,p,4),'-m')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlabel('Days')
% ylabel('u_v (t)')
% xlim([0 10])
% fn = sprintf('Irma_OC_uv_%s',sheet); 
% saveas( gcf, fn, 'png' );
% 
% %Greenville-New Bern-Washington,NC
% p=2;
% figure
% plot(t_f,S_f(p,:,5),'-k+',t_f,S_f(p,:,1),'-ro',t_f,S_f(p,:,2),'-bx',t_f,S_f(p,:,3),'-g^',t_f,S_f(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 5])
% xlabel('Days')
% ylabel('S(t) [%]')
% sheet='GNW';
% fn = sprintf('Irma_OC_S_%s',sheet); 
% saveas( gcf, fn, 'png' );
% 
% figure
% plot(t_f,I_f(p,:,5),'-k+',t_f,I_f(p,:,1),'-ro',t_f,I_f(p,:,2),'-bx',t_f,I_f(p,:,3),'-g^',t_f,I_f(p,:,4),'-md')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% xlim([0 10])
% ylabel('I(t) [%]')
% xlabel('Days')
% fn = sprintf('Irma_OC_I_%s',sheet); 
% saveas( gcf, fn, 'png' )
% 
% 
% figure
% plot(time_f,uvnew_f(:,p,5),'-k',time_f,uvnew_f(:,p,1),'-r',time_f,uvnew_f(:,p,2),'-b',time_f,uvnew_f(:,p,3),'-g',time_f,uvnew_f(:,p,4),'-m')
% legend ('u_v=0','u_v=0.1','u_v=0.25','u_v=0.5','u_v=0.75')
% title('u_{v,max} for Greenville-New Bern-Washington,NC')
% xlabel('Days')
% ylabel('u_v (t)')
% xlim([0 10])

figure
plot(uv1_f,Imax(1,:),'-k+',uv1,Imax(2,:),'-ro',uv1,Imax(3,:),'-bx',uv1,Imax(4,:),'-g^',uv1,Imax(5,:),'-c*',uv1,Imax(6,:),'-md',uv1_f,Imax_f(1,:),'--k',uv1_f,Imax_f(2,:),'-.r')
legend('West Palm' ,'Miami/FT Lauderdale', 'Fort Myers/Naples', 'Tampa/St Pete', 'Orlando', 'Jacksonville','Wilmington (Florence)','Greenville-New Bern-Washington (Florence)')
xlabel('u_{v,max}')
ylabel('Maximum I(t) [%]')
City=['West Palm' ,'Miami/FT Lauderdale', 'Fort Myers/Naples', 'Tampa/St Pete', 'Orlando', 'Jacksonville','Wilmington (Florence)','Greenville-New Bern-Washington (Florence)'];
p=1;
e=[uv1_f',Imax(p,:)]';
xlswrite('Max_Infection_OC.xls', e, City(p), 'A2')
% close all
X=Imax(3,:);
Y=Imax_f(1,:);

uv3=0:0.05:0.35;
uv4=0.35:0.05:1;
y1=-35.528*uv3 + 53.322;
y2=-8.0019*uv4 + 43.723;


uv5=0:0.01:0.39;
uv6=0.39:0.01:1;
y3=-108.26*uv5 + 68.516;
y4=-19.356*uv6 + 33.935;

figure
plot(uv1,Imax_f(1,:),'*r',uv5,y3,'--k',uv6,y4,'--r',0.39*ones(size(0:0.02:26.29,2)),[0:0.02:26.29],'-.r',uv1,Imax(3,:),'or',uv3,y1,'-k',uv4,y2,'--b',0.35*ones(size(0:0.01:40.89,2)),[0:0.01:40.89],'-.b')
legend('Wilmington (Florence)','Bilinear Fcn_1','Bilinear Fcn_2','u_v=0.39','Ft Myers-Naples area','Bilinear Fcn_1','Bilinear Fcn_2','u_v=0,35')
xlabel('u_{v,max}')
ylabel('Maximum I(t) [%]')

