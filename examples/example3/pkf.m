function[]= pkf()
path= './';
file=([path 'hypoDD.reloc']); 
xst= 0.1; yst=0.4; xsh=0.0; ysh=0.20; box=[0.8 0.19]; box1=[0.8 0.14];

icusp=[];

%++++++++++++++++++++++++++++CROSS SECTION++++++++++++++++++++++++++++++++++++
% RELOCATIONS
% ============
axlim= 50;

%=== select for events that are at least mdist from the fault:
mdist= 2.8; 	% dist away from y=m1*x + m2
%--- define fault line:
m(1)= -1.32;	% m1: 1= 45 deg; -1= -45deg...
m(2)= -1.7;	% shift in y-direction

mdat= load(file); 
cusp = mdat(:,1); disp(['events = ' num2str(length(cusp))]);
lon= mdat(:,3); lat=mdat(:,2); 
x0 = mdat(:,5)/1000; y0 = mdat(:,6)/1000; z0 = -mdat(:,4);
ex0 = mdat(:,8)/1000; ey0 = mdat(:,9)/1000; ez0 = mdat(:,10)/1000;
mag=mdat(:,17); 
%ncc= mdat(:,18)+mdat(:,19); nct=mdat(:,20)+mdat(:,21);
yr=mdat(:,11); mo=mdat(:,12); dy=mdat(:,13);
date= yr*10000+mo*100+dy;

%--- get distance of points from fault line:
nv(1)= m(1); nv(2)= m(2); nv(3)= 1; tol= 0.00001;
inc = atan(  sqrt( nv(1)*nv(1) / 1 )  ) *180 / pi; phiA= inc*pi/180;
for i= 1:length(x0);dist(i)= abs(cos(phiA) * ((m(1)*x0(i) + m(2)) - y0(i)));end
ind=find(dist<=mdist);

%--- visually check selections:
figure
plot(x0,y0,'.'); axis('equal'); hold on; plot(x0(ind),y0(ind),'.','color','r');
if(min(x0)<0) yy=  m(1) * 2*min(x0) + m(2); plot([0 2*min(x0)],[m(2) yy]);end
if(max(x0)>0) yy=  m(1) * 2*max(x0) + m(2); plot([0 2*max(x0)],[m(2) yy]);end
hold on; %zoom on; 
figure;
set(gcf,'color',[1 1 1]);
set(gcf,'paperposition',[0.5 0.5 7.5 10],'paperorientation','portrait')

x0=x0(ind);y0=y0(ind);z0=z0(ind);ex0=ex0(ind);ey0=ey0(ind);ez0=ez0(ind); 
mag=mag(ind); cusp=cusp(ind);dist=dist(ind); 
%ncc=ncc(ind);nct=nct(ind);
lon=lon(ind); lat=lat(ind); date= date(ind);

disp(['on fault events = ' num2str(length(cusp))]);

subplot('position',[xst+xsh yst box]); hold on

%=== plot initial locations:
mdat= load([path 'hypoDD.loc']); 
cusp1 = mdat(:,1); 
lon1= mdat(:,3); lat1=mdat(:,2); 
x1 = mdat(:,5)/1000; y1 = mdat(:,6)/1000; z1 = -mdat(:,4);
ex1 = mdat(:,8)/1000; ey1 = mdat(:,9)/1000; ez1 = mdat(:,10)/1000;
mag1=mdat(:,17); 

%--- superimpose relocated events:
for i=1:length(mag);
	nn=20; if(mag(i)<2);nn=6;end
	h= circle(srcrad(mag(i),30)/1000,(x0(i)*cos(phiA)-y0(i)*sin(phiA)),z0(i),'k',nn);
	set(h,'color','k')	
	if(date(i)>20040928); set(h,'color','r');end
end

hold on
text(0,0.6,'MM'); text(17,0.6,'GH');

axis('equal'); axis([-22 35 -16 0]); set(gca,'box','on'); 
xlabel('Distance along strike [km]','fontsize',9); ylabel('Depth [km]','fontsize',9);
set(gca,'yticklabel',[15 10 5 0],'ytick',[-15 -10 -5 0],'fontsize',9)

set(gca,'fontsize',9);


%++++++++++++++++++++++++++++MAP VIEW++++++++++++++++++++++++++++++++++++
subplot('position',[xst yst+ysh box1]); hold on

[th,r]= cart2pol(x1,y1);
th= th+pi/3.5;
[x1,y1]= pol2cart(th,r);
[th,r]= cart2pol(x0,y0);
th= th+pi/3.5;
[x0,y0]= pol2cart(th,r);
for i=1:length(mag);
	h=plot(x0(i),y0(i),'.','markersize',2,'color','k'); hold on
	if(date(i)>=20040928); set(h,'color','r');end
end

axis('equal')

%--- geography and fault trace:
m=load('./sift.lin.1'); a= (m(:,2)-mean(lon))*111.1*cos(mean(lat)*pi/180); b= (m(:,1)-mean(lat))*111.1;
c= m(:,3);
[th,r]= cart2pol(a,b); th= th+pi/3.5; [a,b]= pol2cart(th,r); 
%plot(a,b,'color',[0.6 0.6 0.6],'linewidth',0.3);
plot(a,b,'color','k','linewidth',0.3);
i= find(c==2 |c==3|c==5);
plot(a(i),b(i),'color','k','linewidth',0.3);

%--- stations:
[sta,sla,slo]= textread('hypoDD.sta','%8c %f %f %*[^\n]');
a= (slo-mean(lon))*111.1*cos(mean(lat)*pi/180); b= (sla-mean(lat))*111.1;
[th,r]= cart2pol(a,b); th= th+pi/3.5; [a,b]= pol2cart(th,r); 
plot(a,b,'^','markerfacecolor','k','color','k','markersize',4);

%---SAFOD:
a= (-120.5512-mean(lon))*111.1*cos(mean(lat)*pi/180); b= (35.975-mean(lat))*111.1;
[th,r]= cart2pol(a,b); th= th+pi/3.5; [a,b]= pol2cart(th,r); 
h= plot(a,b,'x','color','k','markersize',12);

axis('equal'); axis([-22 35 -5.9 5.9]); set(gca,'box','on'); 
set(gca,'xticklabel','');
xlabel('','fontsize',9); ylabel('Distance [km]','fontsize',9);
set(gca,'fontsize',9);

title('Parkfield Seismicity, 84-04/10/21; red=04/09/28+')
print -depsc2 Parkfield_84-04.ps

function [rad]= srcrad(magnitude,stressdrop)
%       Matlab function to compute radius of constant stress drop
%       circular source given magnitude and stress drop
%
%       Input:
%               magnitude       earthquake magnitude
%               stress.drop     stress drop in bars
%       Output:
%               radius          source radius in meters

rad = (((7 * 10.^(1.5 * magnitude + 16))./(16 * 1000000 * stressdrop)).^(1/3))/(100);



