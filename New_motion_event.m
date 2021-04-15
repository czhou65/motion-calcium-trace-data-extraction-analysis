function M=New_motion_event(v,lt,ut,flag)%flag=1 motion analysis  0:calcium trace analysis
%{
clear all
close all
[a,b]= uigetfile('*.txt');
cd(b)
kk=strfind(a,'_');
name=a(1:kk(3));
name=strrep(name,'_',' ');
v =Clear_Velocity(a,20);


[lt,ut]=var_low_high_speed(v,20,0.1,2.5,2,2);
%}

%v is the speed trace
%lt/ut is the low/high speed threshold determines by var_low_high_speed.m
%mlt is the locally adjusted low speed threshold
mlt=ones(numel(v),1)*lt;%mlt is initialized with all points equal to lt
ssv=v>lt;%ssv is the initial binary trace with the initial "islands of ones"
nlt=lt;%nlt will be locally adjusted
loop=2;%this process is repeated two times
for k=1:loop
    nlt=nlt+(ut-lt)/10;%increase nlt by 1/10 of the distance between ut and lt    
    A=bwlabel(ssv);%assign indices to each islands
    for i=1:max(A)
        B=find(A==i);%find the location of every island
        if numel(B)>200%only perform adjustment on islands longer than 200 frames
            ssv(B)=v(B)>nlt;%binary trace is locally updated
            mlt(B)=nlt;%local low speed threshold is updated
        end
    end
end
A=bwlabel(ssv);%A contains all the indices and positiion of motion event candidates
MT=1:numel(v);%This motion timestamp is for ploting purpose
M=[];%M will be used to store the onsets and offsets all formal motion events
if flag==1
    ath=(ut-lt)/80;%ath is the acceleration threshold per second per frame
else
    ath=(ut-lt)/60;
end
if flag==1
sl=[20,50];%sl is the search length, in this case, it will search 20 frames prior to 
else
    sl=[20,40];
end
%rough onset and 40 frames after rough onset
%figure(1)
%H1=subplot(2,1,1);
%hold on %these two is for visualization
for i=1:max(A)
    B=find(A==i);% loop through all candidates
    if B(1)~=1%we don't consider event that doesn't start low
        %{
        if flag==1
        M=[M;B(1),B(end)];
        end
        
    else
        %}
        n=numel(B);
        st=mlt(B(1))+(ut-mlt(B(1)))/5;%st is the upperbound speed of the search box
        if B(1)>sl(1) && n>sl(2)%make sure the search box could fit inside the plot
            nB=B(1)-sl(1):B(sl(2)+1);%obtain speed trace indices within the lateral range
            C=diff(v(nB)); %obtain the corresbonding acceleration trace
            if n>30%motion event has to be longer than 30 frames (constrain 3)
                D1=C>ath;%all points that exceed acceleration threshold (constrain 2)
                D2=v(nB)<st;%all points that stay below the upper bound of the search box (constrain 1)
                D=D1&D2(1:end-1);%all points satisfy both constrains
                D3a=bwlabel(D);%loop through all possible region of formal onset
                D3b=bwlabel(D1);%check for consistant acceleration
                D4=zeros(max(D3a),1);
                for j=1:max(D3a)%loop through all possible region of motion onset
                    D5=find(D3a==j);%find indices
                    D6=D3b(D5(1));%indexing back to check for consistant acceleration
                    D4(j)=sum(D3b==D6);%number of points after formal onset that still exceed acceleration threshold
                end
                E=find(D4>=2);%find all region with continuous rapid acceleration of equal to/more than 2 frames
                if ~isempty(E) %if exceed such region
                    for j=1:numel(E)%loop through all formal onset candidate
                        F=find(D3a==E(j));
                        onset=nB(F(1));%onset is the first point of the region
                        C1=diff(v(onset:onset+5));%get the new acceleration trace
                        CC=v(onset:onset+20);
                        C2=diff(CC);
                        C3=find(C2==0);%check for recording artifact
                        C4=CC([C3;C3+1]);
                        if isempty(C4)
                            C4=0;
                        end
                        type=num2str(sum(C1-ath));%represent level of acceleration
                        if sum(C1-ath)>0 && max(C4)<lt%check for constrain 2 and 4
                            if flag==1
                                M=[M;onset,B(end)];
                                break
                            else
                                
                            if onset>41 &&sum(v(onset-41:onset-1)-lt)<10
                            %could be added
                            %to apply stationary constrains prior to onset.
                            M=[M;onset,B(end)];
                            
                            break%formal onset is successfully sellected
                            end
                            end
                        else
                            type='not strong enough/recording bug';
                        end
                    end
                else
                    type='onset too high/not steep enough';
                    onset=B(1);
                end
            else
                onset=B(1);
                type='duration too short'; 
            end
    
      %{
      if onset>23350
        H1=subplot(2,1,1);
        x=-100:100;
        h0=plot(x,v(onset-100:onset+100));
        hold on
        %plot([x(1),x(end)],[lt,lt],'b')
        plot([x(1),x(end)],[ut,ut],'r')
        plot(x,mlt(onset-100:onset+100),'r')
        ylabel('speed (cm/s)')
        %plot([0,0],[0,max(v)],'r')
        %plot([B(end)-onset,B(end)-onset],[0,max(v)],'r')
        pos=[B(1)-onset-sl(1),0,sl(1)+sl(2),st];
        h1=scatter(B(1)-onset,v(B(1)),'r*');
        
        rectangle('Position',pos);
        if B(1)~=onset
            h3=scatter(0,v(onset),'b*');
            legend([h0,h1,h3],{'speed trace','rough onset','tuned onset'})
        else
            legend([h0,h1],{'speed trace','rough onset'})
        end
        hold off

        ylim([0,max(v)])
        xlim([x(1),x(end)])
        title('speed trace near the onset')
        H2=subplot(2,1,2);
        h0=plot(x(1:end-1),diff(v(onset+x(1):onset+x(end)))/20);
        hold on
        h1=plot(B(1)-onset-sl(1):B(1)-onset+sl(2)-1,ones(sl(1)+sl(2),1)*ath/20,'r');
        legend([h0,h1],{'acceleration trace','acceleration threshold'})
        %plot(B(1)-onset-sl(1)+1:B(1)-onset+sl(2)-2,CC,'r')
        hold off
        ylabel('acceleration (cm/s^2)')
        xlabel('frames after tuned onset')
        xlim([x(1),x(end)])
        linkaxes([H1,H2],'x')
        pause  
      end
     %}
    
    %{
     pos=[MT(B(1)-sl(1)),0,sl(1)+sl(2),st];
        scatter(MT(B(1)),v(B(1)),'r*')
        if B(1)~=onset
            scatter(MT(onset),v(onset),'b*')
        end
        rectangle('Position',pos)
        %}
    end
    end
end
%{
plot_motion(MT,v,M);
plot([MT(1),MT(end)],[ut,ut],'r')
plot(MT,mlt,'r')
hold off
ylabel('speed (cm/s)')
title([name,' new method'])
%}
%{
[~,~,M2,~,~]=motion_event(v,lt,ut);
H2=subplot(2,1,2);
plot_motion(MT,v,M2);
plot([MT(1),MT(end)],[ut,ut],'r')
plot([MT(1),MT(end)],[lt,lt],'r')
hold off
xlabel('frames')
title('old method')
linkaxes([H1,H2],'xy')

%}
