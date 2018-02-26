function [out]=coupledlogistic(tslength,r,A,sigma,couplingtype,filename)
%COUPLEDLOGISTIC Generates time-series dynamics for coupled logistic networks of parameter r
%--------------------------------
%Inputs:
%       tslength: length of the time-series (number of points)
%       r: logistic map free parameter. Can be a single number (all nodes
%          are equal) or a vector (specifying different r for each node). 
%       A: adjacency matrix
%       sigma: coupling strength
%       couplingtype: one of the options: 'diffusive' or 'kaneko'
%       filename: name of output file (optional)
%--------------------------------
%Output:
%       out: each column is time-series for a node 
%            described in the adjacency matrix
%--------------------------------
%Usage examples:
%   - Simple X->Y system:
%       A=[0,1;0,0];
%       [x]=coupledlogistic(1e5,4,A,0.2,'diffusive');
%
%
% - Serial:
%       A=[0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 0];
%       [x]=coupledlogistic(1e7,4,A,0.2,'diffusive','serial.mat');
%        
%       (This adjacency matrix A defines a system with nodes [i] connected
%        as: 
%                           [1]->[2]->[3]->[4]. 
%        The dynamics of each node is a r=4 logistic map.
%        Each link is a linear diffusive coupling with strength 0.2.
%        The time-series from each node has 1*10^7 points)
%
% - Parallel:
%       A=[0 1 1 0; 0 0 0 1; 0 0 0 1; 0 0 0 0];
%       [x]=coupledlogistic(1e7,4,A,0.1,'diffusive','parallel.mat');
%        
%       (This adjacency matrix A defines a system with nodes [i] connected
%        as: 
%                               [1]->[2]
%                                |    |
%                                V    V
%                               [3]->[4]. 
%        The dynamics of each node is a r=4 logistic map.
%        Each link is a linear diffusive coupling with strength 0.1.
%        The time-series from each node has 1*10^7 points)
%
% - Wheatstone-bridge:
%       A=[0 1 1 0; 0 0 1 1; 0 0 0 1; 0 0 0 0];
%       [x]=coupledlogistic(1e7,4,A,0.15,'diffusive','wheatstone.mat');
%        
%       (This adjacency matrix A defines a system with nodes [i] connected
%        as: 
%                               [1]->[2]
%                                |  / |
%                                V V  V
%                               [3]->[4]. 
%        The dynamics of each node is a r=4 logistic map.
%        Each link is a linear diffusive coupling with strength 0.15.
%        The time-series from each node has 1*10^7 points)
%--------------------------------
%LaTeX expression:
%
%       If couplingtype='diffusive':
%       $x_{n+1}^i=(1-\sigma)f(x_n^i)+\frac{\sigma}{k_i}\sum_j{A_{ji}(x_n^j-x_n^i)}$
%       Or if couplingtype='kaneko':
%       $x_{n+1}^i=(1-\sigma)f(x_n^i)+\frac{\sigma}{k_i}\sum_j{A_{ji}f(x_n^j)}$
%       where $f(x)=r*x*(1-x)$
%       
%--------------------------------
%(C) Arthur Valencio(1)* and Murilo Baptista(1), 15 December 2017
%   (1)ICSMB, University of Aberdeen,UK
%  *Support: CNPq (Brazil)
    
    disp('Generating time-series');
    nonodes=length(A(1,:));
    
    if length(r)==1
        r(1:nonodes)=r;
    end
    
    
    if strcmp(couplingtype,'diffusive')
        [out]=diffusivecalc(tslength,r,A,sigma,nonodes);
    elseif strcmp(couplingtype,'kaneko')
        [out]=kanekocalc(tslength,r,A,sigma,nonodes);
    end
    
    %cut transient
    out=out(10002:end,:);
    %normalize
    for i=1:nonodes
        out(:,i)=normal(out(:,i));
    end
    %save
    if (exist('filename','var')) 
        disp('saving to file');
        save(filename,'out');
    end
    
end

function [out]=diffusivecalc(tslength,r,A,sigma,nonodes)
%calculation when diffusive
    cond=1;
    count=0;
    while cond    
        temp=0;
        count=count+1;
        %initial cond
        out(1,1:nonodes)=rand(1,nonodes);
        out(2:tslength+10000,1:nonodes)=NaN;
        %calculate
        for n=1:tslength+10000
            %progress bar
            if rem(n,floor(tslength/10))==0
                fprintf('.');
            end
            %actual calc
            for k=1:nonodes
                    %calc coupling    
                    sumterm=0;
                    deg(k)=0;
                    for l=1:nonodes
                        if A(l,k)==1
                            sumterm=sumterm+out(n,l)-out(n,k);
                            deg(k)=deg(k)+1;
                        end
                    end
                    %calc next step
                    if deg(k)>0
                        sumterm=sumterm/deg(k);
                        out(n+1,k)=(1-sigma)*r(k)*out(n,k)*(1-out(n,k))+sigma*sumterm;
                    else %deg=0 means it's an input node, so calc only logistic dynamics
                        out(n+1,k)=r(k)*out(n,k)*(1-out(n,k));
                    end
                    
                    %error: repeat for new initial cond
                    if isnan(out(n+1,k))
                        temp=1;   
                        break;
                    elseif (out(n+1,k)==0)&&(out(n,k)==0)
                        temp=1;
                        break;
                    end                   
            end
        end
        %error handling
        if temp==0
            break;
        else
            disp('recalculating');
        end
        disp(count)
        if count>10
            break;
        end
    end
   
end

function [out]=kanekocalc(tslength,r,A,sigma,nonodes)
%calculation when diffusive
    cond=1;
    count=0;
    while cond
        count=count+1;
        temp=0;
        %initial cond
        out(1,1:nonodes)=rand(1,nonodes);
        out(2:tslength+10000,1:nonodes)=NaN;
        %calculate
        for n=1:tslength+10000
            %progress bar
            if rem(n,floor(tslength/10))==0
                fprintf('.');
            end
            %actual calc
            for k=1:nonodes
                    %calc coupling    
                    sumterm=0;
                    deg(k)=0;
                    for l=1:nonodes
                        if A(l,k)==1
                            sumterm=sumterm+r(k)*out(n,l)*(1-out(n,l));
                            deg(k)=deg(k)+1;
                        end
                    end
                    %calc next step
                    if deg(k)>0
                        sumterm=sumterm/deg(k);
                        out(n+1,k)=(1-sigma)*r(k)*out(n,k)*(1-out(n,k))+sigma*sumterm;
                    else %deg=0 means it's an input node, so calc only logistic dynamics
                        out(n+1,k)=r(k)*out(n,k)*(1-out(n,k));
                    end
                    
                    %error: repeat for new initial cond
                    if isnan(out(n+1,k))
                        temp=1;   
                        break;
                    elseif (out(n+1,k)==0)&&(out(n,k)==0)
                        temp=1;
                        break;
                    end                   
            end
        end
        %error handling
        if temp==0
            break;
        else
            disp('recalculating');
        end
        disp(count)
        if count>10
            break;
        end
    end
    
end

function out=normal(x)
%NORMAL Quick 0 to 1 nomalization
%Arthur Valencio, 6 December 2017

out=(x-min(x))./(max(x)-min(x));

end
