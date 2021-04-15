function [stationary,moving,M,up_transition,down_transition]=motion_event(v,lt,ut)  %608450_01172020
%{
close all
clear all

%% Load plexon and motion files

[AA,b]= uigetfile('*','multiselect','on');
a1=contains(AA,'txt');
a2=contains(AA,'plex');
a=AA{a1};
plex=AA{a2};


cd(b)
%
load(plex)
sti=plx.Stim_onset;
sn=numel(sti);
MT=plx.Timestamp_Motion;
sti(2,:)=plx.Stim_offset;
%

 name=a(1:end-11);

name=strrep(name,'_',' ');

%% Clean up velocity and determine low and high speed thresholds 
%
v =Clear_Velocity(a,20);


[lt,ut]=var_low_high_speed(v,20,0.05,5,2,2);
%
%% determine up transition period


n=numel(v);
MT=MT(1:n);
index=1:n;
%}
n=numel(v);
filter_size=20;
movmean_size=40;
sigmoid_slope=1;
threshold=0.45;
frame_threshold=20;
gap_threshold=20;
filter_trace=zeros(n-filter_size+1,1);
for i=1:numel(filter_trace)
    f=diff(v(i:i+filter_size-1));
    % Why sum(f)?
    filter_trace(i)=sum(f);
end
ft=sigmf(filter_trace,[sigmoid_slope,(ut-lt)/3]);
nft=sigmf(-1*filter_trace,[sigmoid_slope,(ut-lt)/3]);
ft=movmean(ft,movmean_size);
nft=movmean(nft,movmean_size);
% Why 0.45
fft=ft>threshold;
nfft=nft>threshold;
A=bwlabel(fft);
B=bwlabel(nfft);
%{
H1=subplot(2,1,1);
plot(v)
ylabel('speed (cm/s)')
title(['motion trace ',name])
H2=subplot(2,1,2);
plot(ft)
hold on
h1=plot([1,n-20],[threshold,threshold],'r');
legend(h1,'0-1 threshold')
title('filter trace')
xlabel('frames');
ylabel('sigmoid scale')
linkaxes([H1,H2],'x')
%}

%
up_transition=[];
down_transition=[];
for i=1:max(A)
    a=find(A==i);
    t=v(a(1):a(end)+filter_size-1);
    [~,t1]=min(t);
    [~,t2]=max(t);
    %if t2-t1>frame_threshold
        %if isempty(up_transition) || a(1)-up_transition(end,2)>gap_threshold
            % Why a(1)+t1-1,a(1)+t2-1 ?
            up_transition=[up_transition;a(1)+t1-1,a(1)+t2-1];
        %end
    %end
   
end


for i=1:max(B)
    a=find(B==i);
    t=v(a(1):a(end)+filter_size-1);
    [~,t1]=max(t);
    [~,t2]=min(t);
    %if t2-t1>frame_threshold
        %if isempty(up_transition) || a(1)-up_transition(end,2)>gap_threshold
            % Why a(1)+t1-1,a(1)+t2-1 ?
            down_transition=[down_transition;a(1)+t1-1,a(1)+t2-1];
        %end
    %end
   
end
%}

%% Dana's origonal method
a=v>lt;
b=bwlabel(a);
DM=[];
for i=1:max(b)
    c=find(b==i);
    d=v(c);
    e=d>ut;
    f=find(e==1);
    if sum(f)>0
        if f(1)<40
            DM=[DM;c(1),c(end)];
        end
    end
end
    


%% modified Dana's method
SF=0.5; %scaling factor between 0-1  0.5 default, higher, pick up more period
SS=1; %sigmoid slope
%filter_size=20;
sv=sigmf(v,[SS,lt]);
ssv=sv>SF;
A=bwlabel(ssv);

M=[];


for i=1:max(A)
    B=find(A==i);
    %

    %}
    %if B(1)>41 && sum(v(B(1)-41:B(1)-1)-lt)<10
    %
    if B(1)==1
        continue
    
    elseif numel(B)>80 && v(B(41))-v(B(1))>ut-lt
        M=[M;B(1),B(end)];
    elseif numel(B)>40 && v(B(21))-v(B(1))>(ut-lt)/3
        M=[M;B(1),B(end)];
        
        %
    elseif numel(B)>30
        n=numel(B)-10;
        C=zeros(n,1);
        for ii=1:n
            C(ii)=v(B(ii+10))-v(B(ii));
        end
        if max(C)>(ut-lt)/5
            %c1=C>(ut-lt)/6;
            %c2=find(c1==1);
            M=[M;B(1),B(end)];

        end
        %}
    end
    %end
    
        
         

        
   
        
       
        
        

    end


    
    

%% determine stationary periods

stationary=[];
filter_size=40;
stationary_sigmoid_threshold=0.45;
frame_threshold=80;
sigmoid_slope=0.8;


l=1-sigmf(v,[sigmoid_slope,lt]);
L=movmean(l,filter_size);

%{
H1=subplot(2,1,1);
plot(v)
hold on
plot([1,n],[lt,lt],'r')
H2=subplot(2,1,2);
plot(L)
hold on
plot([1,n],[threshold,threshold],'r')
linkaxes([H1,H2], 'x')
%}



L=L>stationary_sigmoid_threshold;
idx=bwlabel(L);
for i=1:max(idx)
       a=find(idx==i);
       N=numel(a);
       if N>frame_threshold 
          stationary=[stationary;a(1),a(end)];
       end
end
%}
%
moving=[];
filter_size=40;
sigmoid_threshold=0.45;
frame_threshold=80;
sigmoid_slope=0.8;


l=sigmf(v,[sigmoid_slope,ut]);
L=movmean(l,filter_size);


%{
H1=subplot(2,1,1);
plot(v)
hold on
plot([1,n],[ut,ut],'r')
H2=subplot(2,1,2);
plot(H)
hold on
plot([1,n],[a2,a2],'r')
linkaxes([H1,H2], 'x')
%}



%
L=L>sigmoid_threshold;
idx=bwlabel(L);
for i=1:max(idx)
       a=find(idx==i);
       N=numel(a);
       if N>frame_threshold 
          moving=[moving;a(1),a(end)];
       end
end
%}
%{
H1=subplot(3,1,1);
[~,~]= plot_motion(MT,v,M);
plot([MT(1),MT(end)],[ut,ut],'r')
hold on
plot([MT(1),MT(end)],[lt,lt],'r')
H2=subplot(3,1,2);
[~,~]= plot_motion(MT,v,m);
H3=subplot(3,1,3);
[~,~]= plot_motion(MT,v,S);
plot([MT(1),MT(end)],[lt,lt],'r')
linkaxes([H1,H2,H3], 'xy')
%}

%{

N=numel(m)/2;
T=randperm(N);
for i=1:10
    subplot(10/2,2,i)
    [~,~]= plot_motion(MT,v,m);
    h3=plot([MT(1),MT(end)],[ut,ut],'r');
    hold on
    plot([MT(1),MT(end)],[lt,lt],'r')
    if i==1
        title(['moving period(2) ',name])
    end
    a1=MT(m(T(i),1));
    a2=MT(m(T(i),2));
    
    xlim([a1-2*(a2-a1),a2+2*(a2-a1)])
end
%xlabel('seconds')
%ylabel('speed (cm/s)')
%legend([h2,h1,h3],{'motion trace','up transition','low/high threshold'})

%}




