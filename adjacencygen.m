function A=adjacencygen(type,nodes)
%ADJACENCYGEN Easy generator of Adjacency Matrix for common directed networks
%--------------------------------
%Inputs:
%       type: 'serial' for serial connection, 'parallel' for parallel,
%             'wheatstone' for wheatstone bridge,'dist3' for distance-3
%       nodes: number of nodes
%--------------------------------
%Output: A: Adjacency Matrix
%--------------------------------
%Example usage:
%       A=adjacencygen('serial',7);
%       is equivalent to [1]->[2]->[3]->[4]->[5]->[6]->[7]
%
%       A=adjacencygen('parallel',4);
%       represents [1]->[2]
%                   |    |
%                   V    V
%                  [3]->[4] 
%
%
%       A=adjacencygen('wheatstone',4);
%       represents [1]->[2]
%                   |  / |
%                   V V  V
%                  [3]->[4]
%
%       A=adjacencygen('dist3',8);
%       represents 
%                  [1]--                   --->[6]
%                       |                  |
%                       |                  |
%                  [2]--->[4]---------->[5]--->[7]
%                       |                  |
%                       |                  |
%                  [3]--                   --->[8]
%--------------------------------
%%(C) Arthur Valencio(1)* and Murilo Baptista(1), 15 December 2017
%(1)ICSMB, University of Aberdeen,UK
%*Support: CNPq, Brazil


if strcmp(type,'serial')
   for i=1:nodes
       for j=1:nodes
           if j==i+1
               A(i,j)=1;
           else
               A(i,j)=0;
           end
       end
   end
   
elseif strcmp(type,'parallel')
    A(1:nodes,1:nodes)=0;
    A(1,2:nodes-1)=1;
    A(2:nodes-1,nodes)=1;
    
elseif strcmp(type,'wheatstone')
    A(1:nodes,1:nodes)=0;
    A(1,2:nodes-1)=1;
    A(2:nodes-1,nodes)=1;
    for i=1:nodes
       for j=1:nodes
           if j==i+1
               A(i,j)=1;
           end
       end
    end

elseif strcmp(type,'dist3')
    A(1:nodes,1:nodes)=0;
    numdendrites=round(nodes/2)-1;
    A(1:numdendrites,numdendrites+1)=1;
    A(numdendrites+1,numdendrites+2)=1;
    A(numdendrites+2,numdendrites+3:nodes)=1;
else
    error('Unknown type');
end


end
