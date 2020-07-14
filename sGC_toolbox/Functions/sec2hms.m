function [TIME,hours, mins, secs] = sec2hms(t)

% Converts time (in seconds) to hours, minutes and seconds.

    hours = floor(t / 3600);
    t = t - hours * 3600;
    mins = floor(t / 60);
    secs = round(t - mins * 60);
    
    space=' ';
    TIMEh=strcat(num2str(hours),'h');
    TIMEm=strcat(num2str(mins),'min');
    TIMEs=strcat(num2str(secs),'secs');
    
    
    if hours>0
        TIME=[TIMEh,space,TIMEm,space,TIMEs];
    elseif mins>0
        TIME=[TIMEm,space,TIMEs];
    else
        TIME=TIMEs;
    end
end