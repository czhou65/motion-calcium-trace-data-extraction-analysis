function v =Clear_Velocity(a,n_moving_average)
%{
close all
clear all
n_moving_average=20;
[a,b]= uigetfile('*.txt');
cd(b)


 name=a(1:end-11);

name=strrep(name,'_',' ');
%}

 f=dir(a);
    Motion_data= importdata(f.name);
    timestamp = Motion_data.data(:,1);
    n=numel(timestamp);
    left_dx = Motion_data.data(:,2);
    left_dy = Motion_data.data(:,3);
    left_dt = Motion_data.data(:,4);
    right_dx = Motion_data.data(:,5);
    water_pin = Motion_data.data(:,8);
    right_dy = Motion_data.data(:,6);
    right_dt = Motion_data.data(:,7);

    ldy=left_dy./left_dt;
    rdy=right_dy./right_dt;
    ldx=left_dx./left_dt;
    rdx=right_dx./right_dt;
    
   
    sensor_angle = 75;
    sr = (sensor_angle/180)*pi;
    
    yl=ldy;
    yr=rdy;
    yl=(yl-yr*cos(sr))/sin(sr);
    v=sqrt(yl.^2+yr.^2)*100;
    v=t_filter(v,[0,100],n_moving_average);   
    
    
   
    
    
    %{
    
    a1=sort(ldy);
    a2=sort(rdy);
    A1=diff(a1)>0.5;
    A2=diff(a2)>0.5;
    A1a=A1(1:round(n/2));
    i1a=find(A1a==1);
    if isempty(i1a)
        lt1=a1(1)-1;
    else
        lt1=a1(i1a(end))+1;
    end
    A1b=[0;A1(round(n/2)+1:end)];
    i1b=find(A1b==1);
    if isempty(i1b)
        ut1=a1(end)+1;
    else
        ut1=a1(i1b(1)+round(n/2))-1;
    end
    Ldy=t_filter(ldy,[lt1,ut1],20);
   
    
    A2a=A2(1:round(n/2));
    i2a=find(A2a==1);
    
     if isempty(i2a)
        lt2=a2(1)-1;
    else
        lt2=a2(i2a(end))+1;
     end
    
    A2b=[0;A2(round(n/2)+1:end)];
    i2b=find(A2b==1);
    if isempty(i2b)
        ut2=a1(end)+1;
    else
         ut2=a2(i2b(1)+round(n/2))-1;
    end
   
    Rdy=t_filter(rdy,[lt2,ut2],20);
    
    
    
    R=(rdx.^2+Rdy.^2).^0.5;
    L=(ldx.^2+Ldy.^2).^0.5;
    V=sqrt(R.^2+L.^2-R.*L*2*cos(sr)/sin(sr)^2)*100;
    V=t_filter(V,[0,100],20);
   
   
    
    %}
 
    
   
    %variance over 2s 40 frames
    %dv=movvar(v,40);
            %{
     H1=subplot(2,1,1);
    plot(1:n,v)
    title('motion trace')
    ylabel('speed cm/s')
    H2=subplot(2,1,2);
    plot(1:n,dv);
    title('moving variance with n=40')
    xlabel('frams')
    linkaxes([H1,H2], 'x')
    %}
    %{
   idx=bwlabel(dv<0.05);
   record=[];
   for i=1:max(idx)
       a=find(idx==i);
       N=numel(a);
       if N>100 && v(a(1))<5
           record=[record;a(1),a(end)];
       end
   end
   record=[record;n,n+1];
   iden=false(n,1);
figure(1)
if numel(record)<3
   plot(1:n,v,'k')
    lt=nan;
    ut=nan;
else
  
  if record(1,1)>1
      I=1:numel(record)/2-1;
      plot(1:record(1,1),v(1:record(1,1)),'k')
      hold on
  else
      plot(record(1,1):record(1,2),v(record(1,1):record(1,2)),'g')
      iden(record(1,1):record(1,2))=true;
      hold on
     plot(record(1,2):record(2,1),v(record(1,2):record(2,1)),'k')
      hold on
      I=2:numel(record)/2-1;
  end
  
  for i=I
       h1=plot(record(i,1):record(i,2),v(record(i,1):record(i,2)),'g');
       iden(record(i,1):record(i,2))=true;
      hold on
      h2=plot(record(i,2):record(i+1,1),v(record(i,2):record(i+1,1)),'k');
      hold on
  end
 
  


vl=v(iden);
vh=v(~iden);
lt=mean(vl)+2*sqrt(var(vl));
ut=mean(vh);

end

h3=plot([1,n],[lt,lt],'r');
hold on
plot([1,n],[ut,ut],'r')
hold off
legend([h1,h2,h3],{'not moving','other','low and high speed threshold'})
title(name)


  
  %}
          
        
  
    
   
    
    

    
    
    
    
   

    
    
   
    
   %{
        
    index=1:n;
    
    
    
    reasonable_vel= distance<100;
    distance=interp1(index(reasonable_vel),distance(reasonable_vel),index);
    distance = movmean(distance,5);
    
    
    
    dtheta=(left_dx+right_dx)/2*100/radius;
    reasonable_vel= (dtheta<pi) & (dtheta>-1*pi);
    dtheta=interp1(index(reasonable_vel),dtheta(reasonable_vel),index);
    dtheta = movmean(dtheta,5);
  
    %}
   

    
    %{
     pathx=zeros(n,1);
    pathy=zeros(n,1);
        dTheta=zeros(n,1);
    for i=2:n
        dTheta(i)=sum(-1*dtheta(1:i));
    end
    
    dY = cos(dTheta).*distance;
    dX = sin(dTheta).*distance;

    for i=2:n
        pathx(i)=sum(dX(1:i));
        pathy(i)=sum(dY(1:i));
    end
    figure(1)
    subplot(4,1,1)
    plot(index,distance)
    title([name,' speed'])
    ylabel('cm/s')
    subplot(4,1,2)
    plot(index,dtheta/pi*180)
    ylim([-180,180])
    ylabel('degrees')
    title('rotation')
    subplot(4,1,[3,4])
    plot(pathx,pathy)
    title('path')
    %}
    %cv=distance';
    %ct=dtheta';
    
   
   
    %{
    
    velocity_cms = sqrt(left_dy.^2+right_dy.^2)*100; 
   
    index  = 1:1:n;
    reasonable_vel= velocity_cms<100;
    velocity_interpolated = interp1(index(reasonable_vel),velocity_cms(reasonable_vel),index);
    clean_velocity = movmean(velocity_interpolated,5);
    cv=clean_velocity';
 
    %}
    