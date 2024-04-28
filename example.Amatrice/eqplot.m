function[]=eqplot(file_reloc,file_sta,file_boot,file_line1,file_line2,phiA,xshift,yshift,box_l,box_w,axlim, minmag, id, itake, ierr, src, year, sta_lab, ev_lab)

if nargin==0;

% PARAMETER SETTING:
path= './';
file_reloc=[path 'hypoDD.reloc']; 
file_loc=[path 'hypoDD.loc']; 
file_sta=[path 'hypoDD.sta']; 
file_boot=[path './bootstrap/fort.999'];
file_boot=[path './bootstrap.out'];
file_line1='/data/hy1/felixw/MATLAB/geo/world.lin';  	% lat, lon (line breaks: NaN)
file_line2='';  				% lat, lon (line breaks: NaN)
file_mrk= ''; 					% lab,lat,lon,elv
file_topo= ''; 					% lat,lon,elv (evenly gridded)

phiA=90; xshift= 0; yshift= 0; box_l= -9; box_w= -9; axlim= -9;
phiA=53; xshift= 0; yshift=  0; box_l=60; box_w=30; axlim= 40;

la1= 30; la2= 34; lo1= 102; lo2= 107;
la1= -9; la2= -9; lo1= -9; lo2= -9;
orlat= -999; orlon= -999;;
orlat= 42.75; orlon= 13.20;;
ierr=0;   		% 0=no errors; 1=SVD errors; 2=bootstrap errors 
			% (from fort.999 or bootstrap.out); 3=both 
src= 0;			% 0= no; e.g. 3= display 3 MPa str dr circ. src 
isymb='.'; symbsiz=1;	% symbol for hypocenters
iorig= 0; iconn= 1;	% plot original locations / connections betw orig-reloc 
ev_lab= 0;		% plot event label
minmag= 1.9;		 
year=0;  		% 1=date in B-B' x-section 
ecc=0;			% ecc> 0: mark events with ecc or more x-corr data 
id=[1344 1801 1958 3278 4828 9825];			% mark these events 
id=[];
itake=[];
sta_lab= 0;		% plot station label
imark= 1; imarklab= 0;	% plot markers from file_mrk
itopo= 0; topofac= 1;	% plot topography (uses file_topo if provided)
idens= 0;               % plot event densitiy distribution >0: box length in km
ifig= [];               % 1= one page, four subplots; 4= four pages.
imap= 1;		% station map projec: 
			% 1=local; 2=regional; 3=global; 0=automatic
end

%========== data processing starts here.... =============================
disp(['Input file: ' char(file_reloc)]);
bx=0;by=0;bz=0;indx=0;
%--- read events 
phiAA= phiA; phiB= phiA-90;  phiB= (phiB-90)*pi/180; phiA= (phiA-90)*pi/180; 

mdat= load(file_reloc); 
cusp = mdat(:,1); mag= mdat(:,17); lon=mdat(:,3); lat=mdat(:,2);
x = mdat(:,5)/1000; y = mdat(:,6)/1000; z = -mdat(:,4);
ex = mdat(:,8)/1000; ey = mdat(:,9)/1000; ez = mdat(:,10)/1000;
yr=mdat(:,11);mo=mdat(:,12);dy=mdat(:,13);
hr=mdat(:,14);mn=mdat(:,15);sc=mdat(:,16);
if(ecc>0); necc= mdat(:,18)+mdat(:,19); else; necc=zeros(1,length(mdat)); end;
if(length(itake) > 0);
	[a,ind,c]= intersect(cusp,itake);
	x= x(ind); y= y(ind); z= z(ind); lon=lon(ind); lat=lat(ind);
	cusp=cusp(ind);mag=mag(ind);ex=ex(ind);ey=ey(ind);ez=ez(ind);
 	yr=yr(ind);mo=mo(ind);dy=dy(ind); hr=hr(ind);mn=mn(ind);sc=sc(ind);
	necc= necc(ind); 
end;

if(orlon==-999); orlon= mean(lon); orlat= mean(lat); end

if(length(unique(cusp))<length(cusp)); disp('Non-unique IDs.'); return;end

%if(abs(sum(diff(x)))<0.001 | length(file_mrk)>0);
	disp('LAT/LON -> KM');
	x=lon2km(orlat,lon-orlon); y=lat2km(lat-orlat);
%end

if(iorig==1); 
	mdat= load(file_loc); 
	lono=mdat(:,3); lato=mdat(:,2);
	xo= mdat(:,5)/1000; yo= mdat(:,6)/1000; zo= -mdat(:,4);cuspo=mdat(:,1);
%	if(abs(sum(diff(xo)))<0.001 | length(file_mrk)>0);
	xo=lon2km(orlat,lono-orlon); yo=lat2km(lato-orlat);
%end
%%	[c,a,b]=intersect(cusp,cuspo);
	[c,b]=ismember(cusp,cuspo);
 	lono=lono(b); lato=lato(b); xo=xo(b); yo=yo(b);zo=zo(b);cuspo=cuspo(b);	
else; 
	xo=0; yo=0; zo=0; cuspo=0; 
end

if(sum(ex)==0); ex=ex+1;ey=ey+1;ez=ez+1;end; mag(find(mag==0))= 0.2;

disp(['Events = ' num2str(length(cusp))]);
disp(['Mean ex = ' num2str(mean(ex))]);
disp(['Mean ey = ' num2str(mean(ey))]);
disp(['Mean ez = ' num2str(mean(ez))]);
disp(['Min mag = ' num2str(min(mag))]);
disp(['Max mag = ' num2str(max(mag))]);
disp(['ORLON = ' num2str(orlon) ';  ORLAT = ' num2str(orlat)]);
if(iorig==1);
	disp(['Mean hypocenter shift = ' num2str(mean(sqrt((x-xo).^2 + (y-yo).^2 + (z-zo).^2))) ' km']);
	disp(['Mean epicenter shift = ' num2str(mean(sqrt((x-xo).^2 + (y-yo).^2))) ' km']);
	disp(['Mean relative hypocenter shift = ' num2str(mean(sqrt(((x - (mean(x)-mean(xo)))-xo).^2 + ((y - (mean(y)-mean(yo)))-yo).^2 + ((z - (mean(z)-mean(zo)))-zo).^2))) ' km']);
	disp(['Mean relative epicenter shift = ' num2str(mean(sqrt(((x - (mean(x)-mean(xo)))-xo).^2 + ((y - (mean(y)-mean(yo)))-yo).^2))) ' km']);
end

%--- Read in station file:
if(length(file_sta>0));
	[sta,slat,slon,dum1,dum2,nccp,nccs,nctp,ncts,rcc,rct,cid]=...
	textread(file_sta,'%8c %f %f %f %f %d %d %d %d %f %f %d');
	ind= find(nctp+ncts+nccp+nccs > 0); 
	slon= slon(ind); slat= slat(ind); sta= sta(ind,:);
	nct=nctp(ind)+ncts(ind); ncc=nccp(ind)+nccs(ind);
else
        sta='';slat=[];slon=[];nccp=[];nccs=[];nctp=[];ncts=[];
	rcc=[];rct=[];cid=[];ncc=[];nct=[];
end
disp(['Stations = ' num2str(length(slat))]);
if(ecc>0); disp(['Events with >= ' num2str(ecc) ' CC data = ' num2str(length(find(necc>=ecc)))]); end

%--- Read line file:
if(length(file_line1)>0);
        geo= load(file_line1);
        xgeo= lon2km(orlat,geo(:,2)-orlon);
        ygeo= lat2km(geo(:,1)-orlat);
        igeo= 1;
else;
        geo=[0,0];
        xgeo= [];
        ygeo= [];
        igeo= 0;
end;
if(length(file_line2)>0);
	faults= load(file_line2); 
	fx= lon2km(orlat,faults(:,2)-orlon);
	fy= lat2km(faults(:,1)-orlat);
	ifault= 1;
else; 
	faults=[0,0];
	fx= []; 
	fy= []; 
	ifault= 0;
end;

%--- Read marker file:
if(length(file_mrk>0));
	[mlab,mla,mlo,mel]=...
	textread(file_mrk,'%s %f %f %f %*[^\n]');
	mx= lon2km(orlat,mlo-orlon);
	my= lat2km(mla-orlat);
else
        mlab='';mla=[];mlo=[];mel=[]; my=[];mx=[]; imark=0;
end

if(la1==-9 | la2==-9 | lo1==-9 |lo2==-9);
        la1a= min([slat' lat']); la2a= max([slat' lat']);
        lo1a= min([slon' lon']); lo2a= max([slon' lon']);
        la1= la1a-(la2a-la1a)/10; la2= la2a+(la2a-la1a)/10;
        lo1= lo1a-(lo2a-lo1a)/10; lo2= lo2a+(lo2a-lo1a)/10;
	if(la1 < -90); la1= -89.9999; end;
	if(la2 >  90); la2=  89.9999; end;
	if(lo1 < -180); lo1= -179.9999; end;
	if(lo2 >  180); lo2=  179.9999; end;

end

%--- Read topo:
if(itopo==1 & length(file_topo)>0);
	topo=load(file_topo);
	k= find(topo(:,1)>lo1 & topo(:,1)<lo2 & topo(:,2)>la1 & topo(:,2)<la2);
	topo= topo(k,:);
	tx= lon2km(orlat,topo(:,1)-orlon);
	ty= lat2km(topo(:,2)-orlat);
	txi= [lon2km(orlat,lo1-orlon):0.1:lon2km(orlat,lo2-orlon)]; 
	tyi= [lat2km(la1-orlat):0.1:lat2km(la2-orlat)]; 
	[xx,yy]=meshgrid(txi,tyi);[xi,yi,zi]=griddata(tx,ty,topo(:,3),xx,yy,'cubic');


        zi=zi/1000;
        tzi=(abs(min(min(zi)))+zi)*topofac;

else
	txi=0; tyi= 0; tzi=0;
end

if(itopo==1 & length(file_topo)==0);
        [latgrat,longrat,Zelv] = satbath(topofac,[la1 la2],[lo1 lo2]);
%       [Z, refvec] =...
%       gtopo30('/home/soho/matlab/topo/gtopo30/',1,[la1 la2],[lo1 lo2]);
%       Z(isnan(Z))= -0.009;
%       [Z1,refvec1]= resizem(Z,4,refvec,'bicubic');
end


%--- Read bootstrap samples:
iboot= 0; ax1= 0; ax2= 0; ax3= 0; az= 0; eid= 0;
if (ierr==2 | ierr==3)
        mdat= load(file_boot); m= size(mdat); 
	if(m(2)==15); iboot= 1; end;
	if(iboot==0);
	        bx= mdat(:,5)/1000; by= mdat(:,6)/1000; bz= mdat(:,7)/1000; indx=mdat(:,1);
%%		fid= fopen('bootstrap.out','w');
%%		fprintf(fid,'YR   MO DY HR MN SEC     LAT      LON       DEPTH   AX1    AX2   AZE-AREA MAG   ID \n');
		ax1= 0; ax2= 0 ; ax3= 0; az= 0; eid= 0;
	else
		ax1= mdat(:,4); ax2= mdat(:,5); ax3= mdat(:,6); az= mdat(:,7); 
		eid= mdat(:,15);
	end
end

if(length(ifig)==0); figure; set(gcf,'clipping','on'); end
if(length(ifig)==1); figure(ifig(1)); set(gcf,'clipping','on'); end

%--- MAP VIEW --- STATION PLOT 
if(length(ifig)==4); figure(ifig(1)); clf; else; subplot(2,2,1); end
hold on

% Local:
%if(length(slat)>0);
if( (imap==0 & (abs(max(slon)-min(slon)) < 10)) | imap==1);
	if(itopo==1 & length(file_topo)>0);
		[xx,yy]=meshgrid(topo(:,1),topo(:,2));
		[xi,yi,zi]=griddata(topo(:,1),topo(:,2),topo(:,3),xx,yy,'cubic');
		pcolor(xi,yi,zi); colormap('jet'); shading flat; hold on;
	end
%	plot(geo(:,1),geo(:,2),'linewidth',0.1,'color',[0 0 0]);
	plot(faults(:,2),faults(:,1),'linewidth',0.1,'color','y');
	if(itopo==1);
        	plot(lon,lat,'.','color','c','markersize',1); else
        	plot(lon,lat,'.','color','r','markersize',1); end
	if(max(ncc)>0); plot(slon(ncc>0),slat(ncc>0),'s','markersize',5,'color','k','markerfacecolor',[1 .7 0],'linewidth',1); end;
	if(max(nct)>0); plot(slon(nct>0),slat(nct>0),'s','markersize',5,'color','k','markerfacecolor','k','linewidth',1); end
	plot(orlon,orlat,'p','markersize',10,'color','k','markerfacecolor','k');
	plot(lon,lat,'.','markersize',3,'color','b','markerfacecolor','b');
	plot(lon(find(mag>=minmag)),lat(find(mag>=minmag)),'o','markersize',10,'color','k');
	if(sta_lab==1);
		for i=1:length(slon); text(slon(i),slat(i),sta(i,:),'fontsize',8);end
	end
	if(imark==1);
        	plot(mlo,mla,'s','color','k','markersize',5);
		if(imarklab==1); text(mlo+0.001,mla,mlab,'color','k','fontsize',9); end
  	end
	axis('equal');box('on');
	axis([lo1 lo2 la1 la2]); grid
	title('Station Map'); xlabel('longitude'); ylabel('latitude'); 

% Regional
elseif((imap==0 & abs(max(slon)-min(slon)) < 40) | imap==2);
%	m=load('coastlines');
	ploc= round((la2-la1)/3*100)/100; mloc= round((lo2-lo1)/3*100)/100;
	axesm('MapProjection','bsam','MapLatLimit',[la1 la2],...
	'MapLonLimit',[lo1 lo2],...
	'MLabelLocation',mloc,'MLineLocation',mloc,'MLabelRound',-2,...
	'PLabelLocation',ploc,'PLineLocation',ploc,'PLabelRound',-2);
	framem; mlabel; plabel; showaxes hide
	if(itopo==1);
            if(length(file_topo)>0);
		[xx,yy]=meshgrid(topo(:,1),topo(:,2));
		[xi,yi,zi]=griddata(topo(:,1),topo(:,2),topo(:,3),xx,yy,'cubic');
		pcolorm(yi,xi,zi); colormap('jet'); shading flat; hold on;
            else
%               Zelv(find(Zelv<0))= NaN;
                surfacem(latgrat,longrat,Zelv);  inc= 1;
                demcmap('inc',Zelv,inc)
            end
	end
	plot3m(faults(:,1),faults(:,2),faults(:,1)*0+11,'color','y');
%	kk= find(m.coastlat*0~=0 | (m.coastlat>min([slat' lat']) &...
%	m.coastlat<max([slat' lat']) & m.coastlon>min([slon' lon']) &...
%	m.coastlon<max([slon' lon'])));
%	plotm(m.coastlat(kk),m.coastlon(kk),'color',[0.5 0.5 0.5]);
	if(max(ncc)>0); plotm(slat(ncc>0),slon(ncc>0),'sk','markersize',5,'markerfacecolor',[1 .7 0],'linewidth',1); end;
	if(max(nct)>0); plotm(slat(nct>0),slon(nct>0),'sk','markersize',5,'linewidth',0.5); end
%	plotm(orlat,orlon,'p','markersize',10,'color','k','markerfacecolor','k');
	plot3m(lat,lon,lon*0+12,'.','markersize',1,'color','b');
	if(length(find(mag>=minmag))>0);
		plot3m(lat(find(mag>=minmag)),lon(find(mag>=minmag)),lat(find(mag>=minmag))*0+12,'o','markersize',10,'color','k');
	end
	if(sta_lab==1);
		for i=1:length(slon); textm(slat(i),slon(i),sta(i,:),'fontsize',8);end
	end
	if(imark==1);
        	plotm(mla,mlo,'s','color','k','markersize',5);
		if(imarklab==1);
        		textm(mla,mlo,char(mlab),'color','k','fontsize',9);
		end
	end

% Global
elseif((imap==0 & abs(max(slon)-min(slon)) >= 40) | imap==3);
        m=load('coast');
        axesm eqaazim; plotm(m.lat,m.long,'color',[0.6 0.6 0.6]);
        setm(gca,'Origin',[orlat orlon]);
        setm(gca,'Flatlimit',180);
        framem; showaxes hide;
	if(itopo==1);
                surfacem(latgrat,longrat,Zelv);  inc= 1;
                demcmap('inc',Zelv,inc)
	end
        plotm(lat,lon,'.','color','r','markersize',1);
	if(max(ncc)>0); plotm(slat(ncc>0),slon(ncc>0),'s','markersize',5,'color','k','markerfacecolor',[1 .7 0],'linewidth',1); end;
	if(max(nct)>0); plotm(slat(nct>0),slon(nct>0),'s','markersize',5,'color','k','markerfacecolor','b','linewidth',1); end
	plotm(orlat,orlon,'p','markersize',10,'color','k','markerfacecolor','k');
	plotm(lat,lon,'.','markersize',2,'color','k','markerfacecolor','k');
        if(sta_lab==1);
                for i=1:length(slon); textm(slat(i),slon(i),sta(i,:),'fontsize',8);end
        end
	if(imark==1);
        	plotm(mla,mlo,'s','color','k','markersize',5);
        	textm(mla,mlo,char(mlab),'color','k','fontsize',9);
	end
end
%end

%--- MAP VIEW --- EVENT PLOT
if(length(ifig)==4); figure(ifig(2)); clf; else; subplot(2,2,2); end
hold on

if(ifault==1); plot(fx,fy,'color',[0.7 0.7 0.7]); end
if(igeo==1); plot(xgeo,ygeo,'color',[0.0 0.0 0.0]); end

if(iorig==1); 
        [i]=ismember(cuspo,cusp);
	plot(xo(i),yo(i),isymb,'markersize',symbsiz,'color',[0.6 0.6 0.6]); 
	if(iconn==1);
		for i=1:length(x);
                        k=find(cusp(i)==cuspo);
			plot([x(i) xo(k)],[y(i) yo(k)],'color',[0.6 0.6 0.6]);
       	 	end
	end
end

plot(x,y,isymb,'markersize',symbsiz,'color','b','linewidth',1);
if(ecc>0); 
	k= find(necc>=ecc);
	plot(x(k),y(k),isymb,'markersize',symbsiz,'color','r','linewidth',1);
end

if(ierr==1 | ierr==3); 
	for i=1:length(x)
		hold on 
		plot([x(i)-ex(i) x(i)+ex(i)],[y(i) y(i)],'color','r');
		plot([x(i) x(i)],[y(i)-ey(i) y(i)+ey(i)],'color','r'); 
	end; 
end

if(ierr==2 | ierr==3);
	if(iboot==0);
		for i=1:length(cusp)
        		ii= find(cusp(i)==indx);
        		if(sum(abs(bx(ii)-mean(bx(ii)))) > 0.00001);
                		[xx,yy,ax1,ax2,az]= conf_ellipse(bx(ii),by(ii),0.9,0);
                		plot(x(i)+xx-mean(bx(ii)),y(i)+yy-mean(by(ii)),'color',[0.7 0.7 0.7]);
%%                		fprintf(fid,'%4i%3i%3i%3i%3i%7.3f%9.4f%10.4f%8.2f%7.2f%7.2f%4.0f%5.1f%7.1f%9i\n',yr(i),mo(i),dy(i),hr(i),mn(i),sc(i),lat(i),lon(i),-z(i),ax1,ax2,az,pi*ax1*ax2,mag(i),cusp(i));
			end
 		end
	else
		for i=1:length(cusp)
			ii= find(cusp(i)==eid);
			[xx,yy]= ellipse(ax1(ii),ax2(ii),az(ii));
                	plot(x(i)+xx,y(i)+yy,'color',[0.7 0.7 0.7]);
		end
	end
end

if(length(id)>0); 
        k= ismember(cusp,id);
        plot(x(k),y(k),'o','markersize',10,'color','g'); 
end;

if(ev_lab==1); 
	for i=1:length(x)
		text(x(i),y(i),num2str(cusp(i)),'fontsize',10,'color','k'); 
	end
end

if(imark==1); plot(mx,my,'s','color','k','markersize',5); end

if(itopo==1);
	%pcolorm(yi,xi,tzi); %%colorbar;
	[c,h]=contour(txi,tyi,tzi); %set(h,'color',[0.7 0.7 0.7]);
end
if(idens>0);
	int= idens;
%	ix=[min(x):(max(x)-min(x))/int:max(x)+(max(x)-min(x))/int];
%	iy=[min(y):((max(y)-min(y))/int):max(y)+(max(y)-min(y))/int];
        ix=[min(x):0.01:max(x)];
        iy=[min(y):0.01:max(y)];

	clear n
	for i=1:length(ix)-1; for j=1:length(iy)-1;
		k=find(x>=ix(i) & x<ix(i+1) & y>=iy(j) & y<iy(j+1));
		n(j,i)= length(k);
%		if(length(k)==0); n(j,i)=NaN; end
	end; end
	pcolor(ix(1:length(ix)-1),iy(1:length(iy)-1),n);
	shading flat
end

plot(x(find(mag>=minmag)),y(find(mag>=minmag)),'o','color','r');
if(axlim==-9); axlim= max(abs([x;y]))*1.3;end
axis('equal'); axis([-axlim axlim -axlim axlim]); set(gca,'box','on');
%title('MAP VIEW'); xlabel('distance [km]'); ylabel('distance [km]');
title(file_reloc); xlabel('distance [km]'); ylabel('distance [km]');
hold on

if(box_l==-9); box_l=axlim*1.8; box_w= axlim*0.9; xshift= 0; yshift= 0; end

%--- plot box location on map 
plot([(-box_l/2)*cos(phiA)-box_w*cos(phiB)+xshift (box_l/2)*cos(phiA)-box_w*cos(phiB)+xshift (box_l/2)*cos(phiA)+box_w*cos(phiB)+xshift (-box_l/2)*cos(phiA)+box_w*cos(phiB)+xshift (-box_l/2)*cos(phiA)-box_w*cos(phiB)+xshift],[(box_l/2)*sin(phiA)+box_w*sin(phiB)+yshift (-box_l/2)*sin(phiA)+box_w*sin(phiB)+yshift (-box_l/2)*sin(phiA)-box_w*sin(phiB)+yshift (box_l/2)*sin(phiA)-box_w*sin(phiB)+yshift (box_l/2)*sin(phiA)+box_w*sin(phiB)+yshift],'color','k');

%--- cross section 
if(length(ifig)==4); figure(ifig(3)); clf; else; subplot(2,2,3); end
hold on

cross_sec(cusp,x,y,z,ex,ey,ez,yr,mo,dy,mag,box_w,box_l,xshift,yshift,ierr,src,year,id,minmag,ev_lab,phiA,phiB,bx,by,bz,indx,iorig,iconn,xo,yo,zo,cuspo,isymb,symbsiz,imark,mx,my,mel,itopo,txi,tyi,tzi,idens,ecc,necc,iboot,ax1,ax2,ax3,az,eid,imap)
set(gca,'box','on'); title('Cross Section: A-A`'); xlabel('distance [km]'); ylabel('depth [km]'); zoom on;

if(length(ifig)==4); figure(ifig(2)); else; subplot(2,2,2); end
hold on
text(-box_l/2*cos(phiA)+xshift, box_l/2*sin(phiA)+yshift,'A'); 
text(box_l/2*cos(phiA)+xshift, -box_l/2*sin(phiA)+yshift,'A`'); 

%--- CROSS SECTION  B-B'
if(length(ifig)==4); figure(ifig(4)); clf; else; subplot(2,2,4); end
hold on
phiA= phiAA+90; tmp= box_w; box_w= box_l/2; box_l= tmp*2; 
phiB= phiA-90;  phiB= (phiB-90)*pi/180; phiA= (phiA-90)*pi/180; 

cross_sec(cusp,x,y,z,ex,ey,ez,yr,mo,dy,mag,box_w,box_l,xshift,yshift,ierr,src,year,id,minmag,ev_lab,phiA,phiB,bx,by,bz,indx,iorig,iconn,xo,yo,zo,cuspo,isymb,symbsiz,imark,mx,my,mel,itopo,txi,tyi,tzi,idens,ecc,necc,iboot,ax1,ax2,ax3,az,eid,imap)
set(gca,'box','on'); title('Cross Section: B-B`'); xlabel('distance [km]'); ylabel('depth [km]'); zoom on;

%--- plot box location on map 
if(length(ifig)==4); figure(ifig(2)); else; subplot(2,2,2); end
hold on
text(-box_l/2*cos(phiA)+xshift, box_l/2*sin(phiA)+yshift,'B'); 
text(box_l/2*cos(phiA)+xshift, -box_l/2*sin(phiA)+yshift,'B`'); 

zoom on
set(gcf,'render','paint');
if(length(ifig)==4);
	figure(ifig(1)); print -depsc eqplot1.ps
	figure(ifig(2)); print -depsc eqplot2.ps
	figure(ifig(3)); print -depsc eqplot3.ps
	figure(ifig(4)); print -depsc eqplot4.ps
else
	print -depsc eqplot.ps
end

%________________________________________________________________________
function[]=cross_sec(cusp,x,y,z,ex,ey,ez,yr,mo,dy,mag,box_w,box_l,xshift,yshift,ierr,src,year,id,minmag,ev_lab,phiA,phiB,bx,by,bz,indx,iorig,iconn,xo,yo,zo,cuspo,isymb,symbsiz,imark,mx,my,mel,itopo,txi,tyi,tzi,idens,ecc,necc,iboot,ax1,ax2,ax3,az,eid,imap)

%i= find(abs((x-xshift)*cos(phiB)-(y-yshift)*sin(phiB))<box_w);
%if(length(i)>0);
%x0=x(i)-xshift;y0=y(i)-yshift;z0=z(i);mag0=mag(i); ex0=ex(i); ey0=ey(i); ez0=ez(i); cusp0=cusp(i); yr0= yr(i); mo0=mo(i); dy0=dy(i); necc0= necc(i);

x0= x-xshift; y0= y-yshift;

i= find(abs((x0)*cos(phiB)-(y0)*sin(phiB)) < box_w &...
abs((x0*cos(phiA)-y0*sin(phiA))) < box_l/2);

if(length(i)>0);
x0=x0(i);y0=y0(i);z0=z(i);mag0=mag(i); ex0=ex(i); ey0=ey(i); ez0=ez(i); 
cusp0=cusp(i); yr0= yr(i); mo0=mo(i); dy0=dy(i); necc0= necc(i);

fid= fopen('events.box','w'); fprintf(fid,'%i\n',cusp0); fclose(fid);

if(itopo==1);
	n= 1;
	for i= 1:length(txi);
		for j= 1:length(tyi);
			txi0(n)= txi(i)-xshift; 
			tyi0(n)= tyi(j)-yshift; 
			tzi0(n)= tzi(j,i); 
			n= n+1;
		end
	end
%	i= find(abs((txi0)*cos(phiB)-(tyi0)*sin(phiB))<box_w);
        i= find(abs((txi0)*cos(phiB)-(tyi0)*sin(phiB)) < box_w &...
        abs((txi0*cos(phiA)-tyi0*sin(phiA))) < box_l/2);
	txi0= txi0(i); tyi0= tyi0(i); tzi0= tzi0(i);
end

if(iorig==1); 
	[i]=ismember(cuspo,cusp0);
	plot3(((xo(i)-xshift)*cos(phiA)-(yo(i)-yshift)*sin(phiA)),zo(i),-ones(length(find(i==1)),1),isymb,'markersize',symbsiz,'color',[0.6 0.6 0.6]); 
	if(iconn==1);
		for i=1:length(x0);
			k=find(cusp0(i)==cuspo);
			plot([(xo(k)-xshift)*cos(phiA)-(yo(k)-yshift)*sin(phiA) x0(i)*cos(phiA)-y0(i)*sin(phiA)],[zo(k) z0(i)],'linewidth',0.1,'color',[0.6 0.6 0.6]); 
        	end
	end
end

if(src==0);
	plot((x0*cos(phiA)-y0*sin(phiA)),z0,isymb,'markersize',symbsiz,'color','b','linewidth',1); 
	if(ecc>0); 
		plot((x0(find(necc0>=ecc))*cos(phiA)-y0(find(necc0>=ecc))*sin(phiA)),z0(find(necc0>=ecc)),isymb,'markersize',symbsiz,'color','r','linewidth',1); 
	end
end

if(ierr== 1 | ierr==3); for i=1:length(x0)
	plot([(x0(i)*cos(phiA)-y0(i)*sin(phiA)-ex0(i)) (x0(i)*cos(phiA)-y0(i)*sin(phiA)+ex0(i))],[z0(i) z0(i)],'color','r');
	plot([(x0(i)*cos(phiA)-y0(i)*sin(phiA)) (x0(i)*cos(phiA)-y0(i)*sin(phiA))],[z0(i)-ez0(i) z0(i)+ez0(i)],'color','r'); 
end;end

if(ierr== 2 | ierr==3);
	if(iboot==0);
		for i=1:length(cusp0)
        		ii= find(cusp0(i)==indx);
        		a=(bx(ii)*cos(phiA)-by(ii)*sin(phiA)); b=bz(ii);
        		[xx,yy,ax1,ax2,az]= conf_ellipse(a,b,0.9,0);
        		plot((x0(i)*cos(phiA)-y0(i)*sin(phiA))-mean(xx)+xx,z0(i)-mean(yy)+yy,'color',[0.7 0.7 0.7])
		end
	else
		for i=1:length(cusp0)
			ii= find(cusp0(i)==eid);
                        [xx,yy]= ellipse(ax1(ii),ax2(ii),az(ii));
                        a= max(xx*cos(phiA)-yy*sin(phiA));
%%        		a= max([(ax1(ii)*cos(phiA)-ax1(ii)*sin(phiA)) (ax2(ii)*cos(phiA)-ax2(ii)*sin(phiA))]);
			[xx,yy]= ellipse(ax3(ii),a,0);
        		plot((x0(i)*cos(phiA)-y0(i)*sin(phiA))-mean(xx)+xx,z0(i)-mean(yy)+yy,'color',[0.7 0.7 0.7])
		end
	end
end

if(src>0); circle(srcrad(mag0,src*10)/1000,x0*cos(phiA)-y0*sin(phiA),z0,'r'); end;

if(year==1); for i=1:length(x0); text((x0(i)*cos(phiA)-y0(i)*sin(phiA)),z0(i),num2str(yr0(i)*10000+mo0(i)*100+dy0(i))); end; end

if(length(id)>0); 
        k= ismember(cusp0,id);
	plot((x0(k)*cos(phiA)-y0(k)*sin(phiA)),z0(k),'o','markersize',10,'color','g'); 
end

if(ev_lab == 1); for i=1:length(x0);hold on
		text((x0(i)*cos(phiA)-y0(i)*sin(phiA)),z0(i),num2str(cusp0(i)),'fontsize',10,'color','k'); end;end

plot((x0(find(mag0>=minmag))*cos(phiA)-y0(find(mag0>=minmag))*sin(phiA)),z0(find(mag0>=minmag)),'o','color','r');

if(imark==1);
	i= find(abs((mx-xshift)*cos(phiB)-(my-yshift)*sin(phiB))<box_w);
	if(length(i)>0);
		plot( ((mx(i)-xshift)*cos(phiA)-(my(i)-yshift)*sin(phiA)),0,'s','markersize',5,'color','k'); end;end;

if(itopo==1 & imap==1);
	tdi= txi0*cos(phiA)-tyi0*sin(phiA);
	d=[min(tdi):(max(tdi)-min(tdi))/20:max(tdi)];
 	for i=1:length(d)-1
		j= find(tdi>=d(i) & tdi<=d(i+1) & tzi0>-9000);
		dd(i)= mean(tzi0(j));
	end
	plot(d(1:length(d)-1)+(d(2)-d(1))/2,dd,'-','color',[0.7 0.7 0.7]);
end

if(idens>0);
	int= idens;
	d=x0*cos(phiA)-y0*sin(phiA);
%	id=[min(d):(max(d)-min(d))/int:max(d)+(max(d)-min(d))/int];
%	iz=[min(z0):((max(z0)-min(z0))/int):max(z0)+(max(z0)-min(z0))/int];
        idst=[min(d):int:max(d)];
        iz=[min(z0):int:max(z0)];

	clear n
	for i=1:length(idst)-1; for j=1:length(iz)-1;
		k=find(d>=idst(i) & d<idst(i+1) & z0>=iz(j) & z0<iz(j+1));
		n(j,i)= length(k);
%		if(length(k)==0); n(j,i)=NaN; end
	end; end
	pcolor(idst(1:length(idst)-1),iz(1:length(iz)-1),n);
	shading flat
end

%axis('equal');axis([-box_l/2 box_l/2 min((min(z0)-(max(z0)-min(z0)+0.01)/5),mean(z0)-box_l/2) max((max(z0)+(max(z0)-min(z0)+0.01)/5),mean(z0)+box_l/2) ]);
axis('equal');axis([-box_l/2 box_l/2 -20 0 ]);

end

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


function h=circle(r,x0,y0,C,Nb)
% CIRCLE adds circles to the current plot
%
% CIRCLE(r,x0,y0) adds a circle of radius r centered at point x0,y0. 
% If r is a vector of length L and x0,y0 scalars, L circles with radii r 
% are added at point x0,y0.
% If r is a scalar and x0,y0 vectors of length M, M circles are with the same 
% radius r are added at the points x0,y0.
% If r, x0,y0 are vector of the same length L=M, M circles are added. (At each
% point one circle).
% if r is a vector of length L and x0,y0 vectors of length M~=L, L*M circles are
% added, at each point x0,y0, L circles of radius r.
%
% CIRCLE(r,x0,y0,C)
% adds circles of color C. C may be a string ('r','b',...) or the RGB value. 
% If no color is specified, it makes automatic use of the colors specified by 
% the axes ColorOrder property. For several circles C may be a vector.
%
% CIRCLE(r,x0,y0,C,Nb), Nb specifies the number of points used to draw the 
% circle. The default value is 300. Nb may be used for each circle individually.
%
% h=CIRCLE(...) returns the handles to the circles.
%
% Try out the following (nice) examples:
%
%% Example 1
%
% clf;
% x=zeros(1,200);
% y=cos(linspace(0,1,200)*4*pi);
% rad=linspace(1,0,200);
% cmap=hot(50);
% circle(rad,x,y,[flipud(cmap);cmap]);
%
%% Example 2
%
% clf;
% the=linspace(0,pi,200); 
% r=cos(5*the);
% circle(0.1,r.*sin(the),r.*cos(the),hsv(40));
% 
%
%% Example 3
%
% clf
% [x,y]=meshdom(1:10,1:10);
% circle([0.5,0.3,0.1],x,y,['r';'y']);
%
%% Example 4
%
% clf
% circle(1:10,0,0,[],3:12);
%
%% Example 5
%
% clf;
% circle((1:10),[0,0,20,20],[0,20,20,0]);

% written by Peter Blattner, Institute of Microtechnology, University of 
% Neuchatel, Switzerland, blattner@imt.unine.ch



% Check the number of input arguments 

if nargin<1,
  r=[];
end;

if nargin==2,
  error('Not enough arguments');
end;

if nargin<3,
  x0=[];
  y0=[];
end;
 

if nargin<4,
  C=[];
end

if nargin<5,
  Nb=[];
end

% set up the default values

if isempty(r),r=1;end;
if isempty(x0),x0=0;end;
if isempty(y0),y0=0;end;
if isempty(Nb),Nb=300;end;
if isempty(C),C=get(gca,'colororder');end;

% work on the variable sizes

x0=x0(:);
y0=y0(:);
r=r(:);
Nb=Nb(:);


if isstr(C),C=C(:);end;

if length(x0)~=length(y0),
  error('length(x0)~=length(y0)');
end;

% how many rings are plottet

if length(r)~=length(x0)
  maxk=length(r)*length(x0);
else
  maxk=length(r);
end;

% drawing loop

for k=1:maxk
  
  if length(x0)==1
    xpos=x0;
    ypos=y0;
    rad=r(k);
  elseif length(r)==1
    xpos=x0(k);
    ypos=y0(k);
    rad=r;
  elseif length(x0)==length(r)
    xpos=x0(k);
    ypos=y0(k);
    rad=r(k);
  else
    rad=r(fix((k-1)/size(x0,1))+1);
    xpos=x0(rem(k-1,size(x0,1))+1);
    ypos=y0(rem(k-1,size(y0,1))+1);
  end;

  the=linspace(0,2*pi,Nb(rem(k-1,size(Nb,1))+1,:)+1);
  h(k)=line(rad*cos(the)+xpos,rad*sin(the)+ypos);
  set(h(k),'color',C(rem(k-1,size(C,1))+1,:));

end;
