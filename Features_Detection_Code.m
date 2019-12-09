%% Data Read and filter
close all;
clear all;
clc;

Data_read_column='B300:B10002'; %N
Data_write_column='7';
peak_distance=00;
data=xlsread('F:\Premier University as Lecturer\Research\Project-2\data\Data Calculation.xlsx',Data_read_column);
data=data;
fs=1000;
n=9;
for i=1:(length(data)-n)
    sum=0;
    for j=i:i+n
       sum=sum+data(j);
    end
    aveg=sum/n;
    data1(i)=aveg; 
end

figure(1)
% plot(data)
% hold on
plot(data1)
grid on

%% r point detection
prompt = 'Put MinPeakHeight: ';
peak_height = input(prompt);
close all;
% [r_point1,r_point_locus1]=findpeaks(data1);
% peak_height=max(r_point1)-0.3;
[r_point,r_point_locus]=findpeaks(data1,'MinPeakHeight',peak_height,'MinPeakDistance',peak_distance);

%%q and s points detection

for i=1:length(r_point)
   for j=r_point_locus(i):-1:2
       diff=data1(j)-data1(j-1);
       if (diff<0)
           q_point(i)=data1(j);
           q_point_locus(i)=j;
           break;
       end
   end
   for j=q_point_locus(i):-1:2
       diff=data1(j)-data1(j-1);
       if (diff>0)
           q_point_start(i)=data1(j);
           q_point_start_locus(i)=j;
           break;
       end
   end
   for j=r_point_locus(i):1:(length(data1)-1)
       diff=data1(j+1)-data1(j);
       if (diff>0)
           s_point(i)=data1(j);
           s_point_locus(i)=j;
           break;
       end
   end
   for j=s_point_locus(i):1:(length(data1)-1)
       diff=data1(j+1)-data1(j);
       if (diff<0)
           s_point_end(i)=data1(j);
           s_point_end_locus(i)=j;
           break;
       end
   end
end

%% t point detection
rr_int_sum=0;
for i=1:(length(r_point_locus)-1)
        rr_int(i)=r_point_locus(i+1)-r_point_locus(i);
        rr_int_sum=rr_int_sum+rr_int(i);
end
rr_int_avg=rr_int_sum/length(rr_int);

for i=1:length(r_point)-1
    x=floor(rr_int_avg/5);
    %x=80;
    X=['subtracted value = ', num2str(x)];
    interval(s_point_end_locus(i):(q_point_start_locus(i+1)-x))=data1(s_point_end_locus(i):(q_point_start_locus(i+1)-x));
    [peaks,locus]=findpeaks(interval);
    max_peak=-100;
    for j=1:length(peaks)
            new_peak=peaks(j);
            if(new_peak>max_peak)
                max_peak=new_peak;
                t_point(i)=max_peak;
                t_point_locus(i)=locus(j);
            end
    end
    interval=0;
end
    disp(X);
    
for i=1:length(t_point)
   for j=t_point_locus(i):-1:2
       diff=data1(j)-data1(j-1);
       if (diff<0)
           t_point_start(i)=data1(j);
           t_point_start_locus(i)=j;
           if(t_point_locus(i)-t_point_start_locus(i)>80)
                break;
           end
       end
   end
   for j=t_point_locus(i):1:length(data1)-1
       diff=data1(j+1)-data1(j);
       if (diff>0)
           t_point_end(i)=data1(j);
           t_point_end_locus(i)=j;
           if(t_point_end_locus(i)-t_point_locus(i)>80)
               break;
           end
       end
   end
end
figure(2)
plot(q_point_locus*(9/9000),q_point,'sr')
hold on
%plot(q_point_start_locus,q_point_start,'*r')
plot(r_point_locus*(9/9000),r_point,'^r')
plot(s_point_locus*(9/9000),s_point,'vr')
%plot(s_point_end_locus,s_point_end,'or')
plot(t_point_locus*(9/9000),t_point,'dr')
%plot(t_point_start_locus,t_point_start,'dr')
%plot(t_point_end_locus,t_point_end,'dr')
plot(((1:length(data1))*(9/9000)),data1)
xlim([0 9])
title('Fetal ECG (arrhythmia)')
xlabel('Time in Seconds','FontSize',10,'FontWeight','bold')
ylabel('Amplitude in mV','FontSize',10,'FontWeight','bold')
grid on
legend('Q point','R point','S point','T point'); 
%% Features

q_point_avg=mean(q_point);
r_point_avg=mean(r_point);
s_point_avg=mean(s_point);
t_point_avg=mean(t_point);

qq_int_sum=0;
for i=1:(length(q_point_locus)-1)
        qq_int(i)=q_point_locus(i+1)-q_point_locus(i);
        qq_int_sum=qq_int_sum+qq_int(i);
end

qq_int_avg=qq_int_sum/length(qq_int);

ss_int_sum=0;
for i=1:(length(s_point_locus)-1)
        ss_int(i)=s_point_locus(i+1)-s_point_locus(i);
        ss_int_sum=ss_int_sum+ss_int(i);
end
ss_int_avg=ss_int_sum/length(ss_int);

qs_int_sum=0;
for i=1:length(s_point_locus)
        qs_int(i)=s_point_locus(i)-q_point_locus(i);
        qs_int_sum=qs_int_sum+qs_int(i);
end
qs_int_avg=qs_int_sum/length(qs_int);

qt_int_sum=0;
for i=1:length(t_point_end_locus)
        qt_int(i)=t_point_end_locus(i)-q_point_start_locus(i);
        qt_int_sum=qt_int_sum+qt_int(i);
end
qt_int_avg=qt_int_sum/length(qt_int);

st_int_sum=0;
for i=1:length(t_point_end_locus)
        st_int(i)=t_point_start_locus(i)-s_point_end_locus(i);
        st_int_sum=st_int_sum+st_int(i);
end
st_int_avg=st_int_sum/length(st_int);

qrs_area_sum=0;
for i=1:length(s_point)
    r_height=r_point(i)-((q_point(i)+s_point(i))/2);
    qrs_area(i)=0.5*r_height*qs_int(i);
    qrs_area_sum=qrs_area_sum+qrs_area(i);
end
qrs_area_avg=qrs_area_sum/length(qrs_area);

energy=0;
for i=1:length(data1)
    energy=energy+(data1(i)*data1(i));
end

average_values=[q_point_avg,r_point_avg,s_point_avg,t_point_avg,qq_int_avg,rr_int_avg,ss_int_avg,qs_int_avg,qt_int_avg,st_int_avg,qrs_area_avg,energy];

qrst=[length(q_point),length(r_point),length(s_point),length(t_point)];
disp(qrst);

xlswrite('F:\Premier University as Lecturer\Research\Project-2\data\Features Data.xlsx',average_values,1,Data_write_column)