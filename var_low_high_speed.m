function [lt,ut,record]=var_low_high_speed(v,Hz,varTH,speedTH,varDur,minDur)
%varTH low threshold for variance to determine stationary motion usually really small  (0.05)
%speedTH  low TH for speed cm/s  (5)
%varDur calculate varient within this time frame in s (2)
%minDur stationary motion has to last for this long in s (2)

a1=varDur*Hz;
a2=minDur*Hz;
a3=varTH;
 n=numel(v);
    
    %variance over 2s a1*Hz frames
  dv=movvar(v,a1);
            
   idx=bwlabel(dv<a3);
   record=[];
   for i=1:max(idx)
       a=find(idx==i);
       N=numel(a);
       % low speed should be less than 5 cm/s and duration should
       % be greater than 2 s - 1- only checking first point in the period
       % v(a(1)<5 - change to mean and check 
       if N>a2 && mean(v(a))<speedTH
           record=[record;a(1),a(end)];
       end
   end
   % Why this n,n+1?
   %record=[record;n,n+1];
   % vector of 0's
   iden=false(n,1);
%figure(1)
% if no stationary period take 10th percentile and if that is > 5 cm/s then
% take 5 cm/s as low speed threshold,
if isempty(record)
   % plot(1:n,v,'k')
    lt=prctile(v,10);
    if lt>2.5
        lt=2.5;
    end
    ut=prctile(v,50);
else
  
% 2-confirm below 


  
  for i=1:numel(record)/2
       iden(record(i,1):record(i,2))=true;
  end
 


vl=v(iden);
vh=v(~iden);
lt=mean(vl)+2*std(vl);
ut=prctile(vh,50);
if ut<5
    ut=5;
end
%{
h3=plot([1,n],[lt,lt],'r');
hold on
plot([1,n],[ut,ut],'r')
hold off
legend([h1,h2,h3],{'not moving','moving','low and high speed threshold'})
%}
end