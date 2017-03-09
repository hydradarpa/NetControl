function [A, pos, A_grid, A_diag] = M_FullTube_prob_cut_add(n,m,p,trials)
% [A, pos] = M_FULLTUBE_PROB_CUT_ADD(n,m,p,trials)
%
%   n:          Number of nodes per slice.
%   m:          Number of slices.
%   p:          Unknown. Defaults to (1/m)*ones(m,1) to generate uniform probabilities.
%   trials:     Unknown. Defaults 10.
%
% Returns
%   A:          Adjacency matrix of size [n*m, n*m]
%   pos:        3D positions of each node 1,...,n*m, where each column is 1 position.
%   A_grid:     Portion of 'A' that corresponds to the 'grid' edges.
%   A_diag:     Portion of 'A' that corresponds to the diagonal cross edges.

% CREATES the adjacency matrix of a size n by m, where n - rows, m -
% columns

% CUTS the connections to fit the probability distrubution p.
% ADD the connections to restore the impaired the probability distribution
% p - vector of the outcomes consisting of probabilites

if (nargin < 3); %p = (1/m)*ones(m,1);
p = [0, 0.0714, 0.6429, 0.2857]; end;
if (nargin < 4); trials = 100; end;


function ind = graph_element(x,y,m)
%UNTITLED2 Summary of this function goes here
%   returns the number of an element on the graph, numbrering: 
% left -> right, up -> down
ind=x+(y-1)*m;
end

% N - size parameter, number of elements in the matrix - sqrt(elements)


A_grid = zeros(m*n,m*n);
A_diag = zeros(m*n,m*n);

%
% Generate node posisions
%
pos = zeros(3, m*n);
dTheta = (2*pi)/n;
r = dTheta;
idx = 1;
for iTheta = 0:(n-1)
    x = r*cos(iTheta*dTheta);
    y = r*sin(iTheta*dTheta);
        
    for iZ = 0:(m-1)
        z = iZ*r;
        
        pos(:,idx) = [x;y;z];
        idx = idx+1;
    end
end

for x=1:1:m         % loop over all all graph edges, generate tube-like adjacency matrix
    for y=1:1:n
                         
    % tube border conditions    
        A_grid(graph_element(x,1,m),graph_element(x,n,m))=1;        
        A_grid(graph_element(x,n,m),graph_element(x,1,m))=1;                
                               
        
        if x>1                                                     % diagonal border
            A_diag(graph_element(x,1,m),graph_element(x-1,n,m))=1;        
            A_diag(graph_element(x-1,n,m),graph_element(x,1,m))=1;                
        end
        
        if x+1<=m                                                  % diagonal border
            A_diag(graph_element(x,1,m),graph_element(x+1,n,m))=1;        
            A_diag(graph_element(x+1,n,m),graph_element(x,1,m))=1;                
        end
        
        A_diag(graph_element(1,y,m),graph_element(m,y,m))=0;                        
        A_diag(graph_element(m,y,m),graph_element(1,y,m))=0;        
                    
    % formula for the rest of the elements
        
       if x>1
       A_grid(graph_element(x,y,m),graph_element(x-1,y,m))=1;            
       A_grid(graph_element(x-1,y,m),graph_element(x,y,m))=1;
        if y>1                                                      % diagonal elements
            A_diag(graph_element(x,y,m),graph_element(x-1,y-1,m))=1;              
            A_diag(graph_element(x-1,y-1,m),graph_element(x,y,m))=1;
        end
        if y+1<=n                                                   % diagonal elements
            A_diag(graph_element(x,y,m),graph_element(x-1,y+1,m))=1;            
            A_diag(graph_element(x-1,y+1,m),graph_element(x,y,m))=1;
        end
       end       
   
       
       if y>1      
        A_grid(graph_element(x,y,m),graph_element(x,y-1,m))=1;
        A_grid(graph_element(x,y-1,m),graph_element(x,y,m))=1;        
       end                                
           
       
       if x+1<=m
       A_grid(graph_element(x,y,m),graph_element(x+1,y,m))=1;                   
       A_grid(graph_element(x+1,y,m),graph_element(x,y,m))=1;
           if y>1                                                   % diagonal elements
            A_diag(graph_element(x,y,m),graph_element(x+1,y-1,m))=1;            
            A_diag(graph_element(x+1,y-1,m),graph_element(x,y,m))=1;
           end                                                     
            if y+1<=n                                               % diagonal elements
            A_diag(graph_element(x,y,m),graph_element(x+1,y+1,m))=1;            
            A_diag(graph_element(x+1,y+1,m),graph_element(x,y,m))=1;
            end
       end
      
       if y+1<=n
       A_grid(graph_element(x,y,m),graph_element(x,y+1,m))=1;
       A_grid(graph_element(x,y+1,m),graph_element(x,y,m))=1;              
       end
                            
       
    end
end

A = A_grid + A_diag;

% counter
q=0;


% Randomly cutting connections

for k=1:1:100           % amount of trials to find full connected graph

p_connect = makedist('Multinomial','Probabilities',p);
B=A;                        % replace the connectivity


p_rand=random(p_connect,1,m*n); % create vector of connection number per element
     
for kk=1:1:trials

% Remove connections

for i=1:1:m*n                   % loop over all all element, half of the symmetric matrix
    
    connected=find(B(i,:)>0);   % indexes of connected elements        
    % remove connections
    if length(connected) >= p_rand(i) % cut elements if there too many connections        
    connection_cut=datasample(connected,abs(length(connected)-p_rand(i)),'Replace',false);         
    B(i,connection_cut)=0;
    B(connection_cut,i)=0;
    end
    
end
 
% Add missing connections


% add more elements if too many are cut    
for i=1:1:m*n                   % loop over all the elements, half of the symmetric matrix
    
    connected=find(B(i,:)>0);   % indexes of connected elements
        
    if length(connected) < p_rand(i) % add elements if there are not enough connections
    connection_possible=find(A(i,:)>0);    
    connection_add=datasample(connection_possible,abs(length(connected)-p_rand(i)),'Replace',false);
    B(i,connection_add)=1;
    B(connection_add,i)=1;
    end
    
end

for i=1:1:m*n                   % loop over all all element, half of the symmetric matrix
    
    connected=find(B(i,:)>0);   % indexes of connected elements        
    % remove connections
    if length(connected) >= p_rand(i) % cut elements if there too many connections        
    connection_cut=datasample(connected,abs(length(connected)-p_rand(i)),'Replace',false);         
    B(i,connection_cut)=0;
    B(connection_cut,i)=0;
    end
    
end


end

%}

B_graph=graph(B);       % create the graph out of the adjacency matrix


if max(conncomp(B_graph))==1    
    break              % break the loop if all elements are parts of the same graph
else
    disp('Did not find full connected graph yet');
end


end

A=B;

end
