function[]=eqplot
% Matlab script to read and display hypoDD output files.
% Edit parameters at the beginning of the script to make changes
% to the view of the plot.

% PARAMETER DESCRIPTION:
%--- file names:
% file_loc      hypocenter file (hypoDD.loc or hypoDD.reloc)
% file_sta      station file (hypoDD.sta)
% file_line     line file to display geography
%                 (lines of lon, lat; separators are Nan, NaN)

%--- cross section:
% phiA          strike of cross section
% xshift        x-location
% yshift        y-location
% box_l         length
% box_w         half width

%--- event selection:
% axlim         axis dimension of map view plot
% minmag        circle events with mag>minmag
% id=[]         mark these events
% itake=[]      display only these events
% err=1         display errors [0,1]

% PARAMETER SETTING: 
file_loc='./hypoDD.reloc'; file_sta='hypoDD.sta'; file_line='';
phiA=140; xshift=0; yshift=0; box_l=14; box_w=3;
axlim=7; minmag=10; id=[]; itake=[]; err=0;
 

%========== data processing starts here....
%--- read events 
phiAA= phiA; phiB= phiA-90;  phiB= (phiB-90)*pi/180; phiA= (phiA-90)*pi/180; 
mdat=load(file_loc); 
cusp = mdat(:,1); mag= mdat(:,17); lon=mdat(:,3); lat=mdat(:,2);
x = mdat(:,5)/1000; y = mdat(:,6)/1000; z = -mdat(:,4);
ex = mdat(:,8)/1000; ey = mdat(:,9)/1000; ez = mdat(:,10)/1000;
if(sum(ex)==0); ex=ex+1;ey=ey+1;ez=ez+1;end; mag(find(mag==0))= 0.2;
if(length(itake) > 0);
	clear ind; k=1; for i=1:length(cusp); for j=1:length(itake); 
	if(itake(j)==cusp(i)); ind(k)= i; k=k+1; end;end;end;
	x= x(ind); y= y(ind); z= z(ind); lon=lon(ind); lat=lat(ind);
	cusp=cusp(ind);mag=mag(ind);ex=ex(ind);ey=ey(ind);ez=ez(ind);
end;
disp(['# of events = ' num2str(length(cusp))]);
disp(['mean ex = ' num2str(mean(ex))]);
disp(['mean ey = ' num2str(mean(ey))]);
disp(['mean ez = ' num2str(mean(ez))]);
disp(['min mag = ' num2str(min(mag))]);
disp(['max mag = ' num2str(max(mag))]);
%length(cusp(find(z0>-32.3)))'

%--- read stations:
fid1=fopen(file_sta,'r');
for i=1:1000000
        [dum,count]= fscanf(fid1,'%s1/');
        if count == 0 break; end
        sta(i,1:length(dum))= dum;
        slat(i)=fscanf(fid1,'%f1') ; slon(i)=fscanf(fid1,'%f1') ;
        dum=fscanf(fid1,'%f1') ; dum=fscanf(fid1,'%f1') ;
        np(i)=fscanf(fid1,'%d1') ; ns(i)=fscanf(fid1,'%d1') ;
        nnp(i)=fscanf(fid1,'%d1') ; nns(i)=fscanf(fid1,'%d1') ;
        dum=fscanf(fid1,'%f1'); dum=fscanf(fid1,'%f1'); dum=fscanf(fid1,'%d1/');
end
ind= find(np+ns+nnp+nns > 0); slon= slon(ind); slat= slat(ind); sta= sta(ind,:);
disp(['# of stations = ' num2str(length(sta))]);

%--- read lines:
if(length(file_line)>0);border=load(file_line); else; border=[0,0];end;

figure; set(gcf,'clipping','on');

%--- STATION PLOT  (map view)
subplot(2,2,1);hold on; 
plot(border(:,1),border(:,2),'linewidth',0.1,'color',[0.8 0.8 0.8]);
plot(lon,lat,'o','markersize',2,'color','r');
plot(slon,slat,'v','markersize',2); 
for i=1:length(slon); text(slon(i),slat(i),sta(i,:),'fontsize',8);end
axis([min([slon lon']) max([slon lon']) min([slat lat']) max([slat lat'])]); 
title('STATION MAP'); xlabel('longitude'); ylabel('latitude'); 
axis('equal');box('on');

%--- PLOT EVENTS (map view)
subplot(2,2,2);hold on
plot(x,y,'.','markersize',1,'color','r');
if(err==1); for i=1:length(x)
	hold on 
	plot([x(i)-ex(i) x(i)+ex(i)],[y(i) y(i)],'color','r');
	plot([x(i) x(i)],[y(i)-ey(i) y(i)+ey(i)],'color','r'); end; end
if(length(id)>0); for i=1:length(x)
	hold on 
	for k=1:length(id); if(id(k)==cusp(i)); 
		plot(x(i),y(i),'o','markersize',10,'color','g'); end;end;
end;end;
plot(x(find(mag>minmag)),y(find(mag>minmag)),'o','color','r');
axis('equal'); axis([-axlim axlim -axlim axlim]); set(gca,'box','on');
title('MAP VIEW'); xlabel('distance [km]'); ylabel('distance [km]');
hold on

%--- CROSS SECTION  A-A'
%--- plot box location on map 
subplot(2,2,2); hold on
plot([(-box_l/2)*cos(phiA)-box_w*cos(phiB)+xshift (box_l/2)*cos(phiA)-box_w*cos(phiB)+xshift (box_l/2)*cos(phiA)+box_w*cos(phiB)+xshift (-box_l/2)*cos(phiA)+box_w*cos(phiB)+xshift (-box_l/2)*cos(phiA)-box_w*cos(phiB)+xshift],[(box_l/2)*sin(phiA)+box_w*sin(phiB)+yshift (-box_l/2)*sin(phiA)+box_w*sin(phiB)+yshift (-box_l/2)*sin(phiA)-box_w*sin(phiB)+yshift (box_l/2)*sin(phiA)-box_w*sin(phiB)+yshift (box_l/2)*sin(phiA)+box_w*sin(phiB)+yshift]);
text(-box_l/2*cos(phiA)+xshift, box_l/2*sin(phiA)+yshift,'A'); 
text(box_l/2*cos(phiA)+xshift, -box_l/2*sin(phiA)+yshift,'A`'); 

%--- cross section 
subplot(2,2,3); hold on
i= find(abs((x-xshift)*cos(phiB)-(y-yshift)*sin(phiB))<box_w);
if(length(i)>0);
x0=x(i)-xshift;y0=y(i)-yshift;z0=z(i);mag0=mag(i); ex0=ex(i); ey0=ey(i); ez0=ez(i); cusp0=cusp(i);
plot((x0*cos(phiA)-y0*sin(phiA)),z0,'.','markersize',1,'color','r');
if(err== 1); for i=1:length(x0)
	plot([(x0(i)*cos(phiA)-y0(i)*sin(phiA)-ex0(i)) (x0(i)*cos(phiA)-y0(i)*sin(phiA)+ex0(i))],[z0(i) z0(i)],'color','r');
	plot([(x0(i)*cos(phiA)-y0(i)*sin(phiA)) (x0(i)*cos(phiA)-y0(i)*sin(phiA))],[z0(i)-ez0(i) z0(i)+ez0(i)],'color','r'); end;end
if(length(id)>0); for i=1:length(x0);hold on
	for k=1:length(id); if(id(k)==cusp0(i)); 
		plot((x0(i)*cos(phiA)-y0(i)*sin(phiA)),z0(i),'o','markersize',10,'color','g'); end;end
end;end
plot((x0(find(mag0>minmag))*cos(phiA)-y0(find(mag0>minmag))*sin(phiA)),z0(find(mag0>minmag)),'o','color','r');
axis('equal');axis([-box_l/2 box_l/2 min((min(z0)-(max(z0)-min(z0)+0.01)/5),mean(z0)-box_l/2) max((max(z0)+(max(z0)-min(z0)+0.01)/5),mean(z0)+box_l/2) ]);
set(gca,'box','on'); title('Cross Section: A-A`'); xlabel('distance [km]'); ylabel('depth [km]'); %zoom on;
end;

%--- CROSS SECTION  B-B'
phiA= phiAA+90; tmp= box_w; box_w= box_l/2; box_l= tmp*2; 
phiB= phiA-90;  phiB= (phiB-90)*pi/180; phiA= (phiA-90)*pi/180; 

%--- plot box location on map 
subplot(2,2,2); hold on
plot([(-box_l/2)*cos(phiA)-box_w*cos(phiB)+xshift (box_l/2)*cos(phiA)-box_w*cos(phiB)+xshift (box_l/2)*cos(phiA)+box_w*cos(phiB)+xshift (-box_l/2)*cos(phiA)+box_w*cos(phiB)+xshift (-box_l/2)*cos(phiA)-box_w*cos(phiB)+xshift],[(box_l/2)*sin(phiA)+box_w*sin(phiB)+yshift (-box_l/2)*sin(phiA)+box_w*sin(phiB)+yshift (-box_l/2)*sin(phiA)-box_w*sin(phiB)+yshift (box_l/2)*sin(phiA)-box_w*sin(phiB)+yshift (box_l/2)*sin(phiA)+box_w*sin(phiB)+yshift]);
text(-box_l/2*cos(phiA)+xshift, box_l/2*sin(phiA)+yshift,'B'); 
text(box_l/2*cos(phiA)+xshift, -box_l/2*sin(phiA)+yshift,'B`'); 

%--- cross section 
subplot(2,2,4); hold on
i= find(abs((x-xshift)*cos(phiB)-(y-yshift)*sin(phiB))<box_w);
if(length(i)>0);
x0=x(i)-xshift;y0=y(i)-yshift;z0=z(i);mag0=mag(i); ex0=ex(i); ey0=ey(i); ez0=ez(i); cusp0=cusp(i);
plot((x0*cos(phiA)-y0*sin(phiA)),z0,'.','markersize',1,'color','r');
if(err== 1); for i=1:length(x0)
	plot([(x0(i)*cos(phiA)-y0(i)*sin(phiA)-ex0(i)) (x0(i)*cos(phiA)-y0(i)*sin(phiA)+ex0(i))],[z0(i) z0(i)],'color','r');
	plot([(x0(i)*cos(phiA)-y0(i)*sin(phiA)) (x0(i)*cos(phiA)-y0(i)*sin(phiA))],[z0(i)-ez0(i) z0(i)+ez0(i)],'color','r'); end;end
if(length(id)>0); for i=1:length(x0);hold on
	for k=1:length(id); if(id(k)==cusp0(i)); 
		plot((x0(i)*cos(phiA)-y0(i)*sin(phiA)),z0(i),'o','markersize',10,'color','g'); end;end
end;end
plot((x0(find(mag0>minmag))*cos(phiA)-y0(find(mag0>minmag))*sin(phiA)),z0(find(mag0>minmag)),'o','color','r');
axis('equal');axis([-box_l/2 box_l/2 min((min(z0)-(max(z0)-min(z0)+0.01)/5),mean(z0)-box_l/2) max((max(z0)+(max(z0)-min(z0)+0.01)/5),mean(z0)+box_l/2) ]);
set(gca,'box','on'); title('Cross Section: B-B`'); xlabel('distance [km]'); ylabel('depth [km]'); %zoom on;
end;

