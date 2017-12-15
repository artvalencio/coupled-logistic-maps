clear all;
limnodes=25;
coupling=0.3;
npoints=10000;

figure(1)
for i=3:limnodes
    disp(i)
    A=adjacencygen('serial',i);
    [x,L,pL]=coupledlogistic(npoints,4,A,coupling,'diffusive',strcat('serial',num2str(i),'.mat'));
    cor(i)=correlation(x(:,1),x(:,end));
    lap(i)=pL(1,end);
end
scatter(lap(3:limnodes),cor(3:limnodes),25,3:limnodes);
xlabel('Value of pseudo-inverse of the Laplacian between first and last node');
ylabel('Correlation between first and last node');
title('Growing serial network')

figure(2)
for i=3:limnodes
    cor(i)=correlation(x(:,1),x(:,i));
end
scatter(pL(1,3:limnodes),cor(3:limnodes),25,3:limnodes)
xlabel('Value of pseudo-inverse of the Laplacian between first and nth node')
title(strcat(num2str(limnodes),'node serial network'));
ylabel('Correlation between first and nth node')


figure(3)
for i=4:limnodes
    disp(i)
    A=adjacencygen('parallel',i);
    [x,L,pL]=coupledlogistic(npoints,4,A,coupling,'diffusive',strcat('parallel',num2str(i),'.mat'));
    cor(i)=correlation(x(:,1),x(:,end));
    lap(i)=pL(1,end);
end
scatter(lap(4:limnodes),cor(4:limnodes),25,4:limnodes);
xlabel('Value of pseudo-inverse of the Laplacian between first and last node');
ylabel('Correlation between first and last node');
title('Growing parallel network')

figure(4)
for i=4:limnodes
    cor(i)=correlation(x(:,1),x(:,i));
end
scatter(pL(1,4:limnodes),cor(4:limnodes),25,4:limnodes)
xlabel('Value of pseudo-inverse of the Laplacian between first and nth node')
title(strcat(num2str(limnodes),'node parallel network'))
ylabel('Correlation between first and nth node')

figure(5)
for i=4:limnodes
    disp(i)
    A=adjacencygen('wheatstone',i);
    [x,L,pL]=coupledlogistic(npoints,4,A,coupling,'diffusive',strcat('wheatstone',num2str(i),'.mat'));
    cor(i)=correlation(x(:,1),x(:,end));
    lap(i)=pL(1,end);
end
scatter(lap(4:limnodes),cor(4:limnodes),25,4:limnodes);
xlabel('Value of pseudo-inverse of the Laplacian between first and last node');
ylabel('Correlation between first and last node');
title('Growing wheatstone network')

figure(6)
for i=4:limnodes
    cor(i)=correlation(x(:,1),x(:,i));
end
scatter(pL(1,4:limnodes),cor(4:limnodes),25,4:limnodes)
xlabel('Value of pseudo-inverse of the Laplacian between first and nth node')
title(strcat(num2str(limnodes),'node wheatstone network'));
ylabel('Correlation between first and nth node')
