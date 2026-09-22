% This function is used to find the potential trend turning point in long
% term datsets.
% ***************************************************************
% Code written by Pankaj Dey, PhD scholar, Department of Civil Engineering,
% Indian Institute of Science, Bangalore.
% ***************************************************************
%
%
%  Input: The time series data with first column with time index and the
%  second column with the values of the variable. Size is n x 2, where n is
%  no. of elements in time series dataset.
%
%
%  Output: 'sqmk' consists of 3 columns. The first column gives the
%  progressive series, second column gives the retrograde series and the
%  third one gives the logical series. The ones are the value in the time
%  index from where the potential truning point with respect ot trend has
%  occured.
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Reference: Chatterjee, S., Bisai, D. and Khan, A., 2014. Detection of approximate potential trend turning points in temperature time series (1941�2010) for Asansol Weather Observation Station, West Bengal, India. India. Atmospheric and Climate Sciences, 4, pp.64-69.
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sqmk]=seqMK(data)
% Progressive series computation
year=data(:,1);
m=length(data);
for j=1:m
    for i=1:j
        if (data(j,2)>data(i,2))
            n(j,i)=1;
        else
            n(j,i)=0;
    end
    end
end
for i=1:m
    num(i,1)=sum(n(i,:));
    t=cumsum(num);
    var(i,1)=i*(i-1)*(2*i+5)/72;
    mean(i,1)=i*(i-1)/4;
    u_prog(i,1)=(t(i,1)-mean(i,1))/(sqrt(var(i,1)));  % it is the progressive series
end
% Retrograde series computation
data1=flipud(data(:,2));  % inverting the entire dataset
for j=1:m
    for i=1:j
        if (data1(j,1)>data1(i,1))
            n1(j,i)=1;
        else
            n1(j,i)=0;
    end
    end
end
for i=1:m
    num1(i,1)=sum(n1(i,:));
    t1=cumsum(num1);
    var1(i,1)=i*(i-1)*(2*i+5)/72;
    mean1(i,1)=i*(i-1)/4;
    u_retr(i,1)=(t1(i,1)-mean1(i,1))/(sqrt(var1(i,1)));
end
u_retr=flipud(u_retr); % inverting the retrograde series
% Estimation of the turning points
diff=u_prog-u_retr;
for i=2:length(diff)-2
    if(sign(diff(i))==sign(diff(i+1)))
        turning_point(i+1,1)=0;
    else
        turning_point(i+1,1)=1;
    end
end
turning_point(m,1)=0;
% Final output array consists of the progressive series, retrograde series
% and the turning points.
sqmk=[u_prog, u_retr,turning_point];
% Plotting the progressive and the retrograde series
% plot(year,u_prog,'r')
% hold on
% plot(year,u_retr,'b')
% grid on
% title('Sequential Mann Kendall Test');
% xlabel('Time');
% legend( 'Progressive Series', 'Retrograde Series');