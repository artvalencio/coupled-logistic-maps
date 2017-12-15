function corr = correlation(x,y)
%Calculates correlation the Nicolas way

    corr=0;
    avg1=mean(x);
    avg2=mean(y);
    for i=1:length(x)
            corr=corr+(x(i)-avg1).*(y(i)-avg2);
    end
    corr=corr/length(x);

end