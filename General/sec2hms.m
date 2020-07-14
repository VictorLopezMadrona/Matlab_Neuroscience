function [TIME,hours, mins, secs] = sec2hms(t)
    hours = floor(t / 3600);
    t = t - hours * 3600;
    mins = floor(t / 60);
    secs = round(t - mins * 60);
    
    espacio=' ';
    TIMEh=strcat(num2str(hours),'h');
    TIMEm=strcat(num2str(mins),'min');
    TIMEs=strcat(num2str(secs),'secs');
    
    
    if hours>0
        TIME=[TIMEh,espacio,TIMEm,espacio,TIMEs];
    elseif mins>0
        TIME=[TIMEm,espacio,TIMEs];
    else
        TIME=TIMEs;
    end
end