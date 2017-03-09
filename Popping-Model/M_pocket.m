function [A,pos] = M_pocket(n,m,p0,sigma_x,sigma_y)
% generates the spatially dependent matrix on the tube graph
% euclidean distance. Probability distribution, sigma_x, sigma_y

% output A - one instance of generated matrix, A_prob - prob matrix

if (nargin < 3); p0 = 1; end;
if (nargin < 4); sigma_x = 1.5; end;
if (nargin < 5); sigma_y = 1.5; end;

function ind = graph_element(x,y,m)
%   returns the number of an element on the graph, numbrering: 
% left -> right, up -> down
ind=x+(y-1)*m;
end


A=zeros(m*n,m*n);
A_prob=zeros(m*n,m*n);

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

% Now these connections are just a sheet of paper
% Need to add border conditions

for x=1:1:m         % loop over all nodes
    for y=1:1:n

        for i=1:1:m % loop over all connections
            for j=1:1:n
                
                % internal elements
                p_random=rand;
                p=p0*exp(-(x-i)^2/sigma_x -(y-j)^2/sigma_y);    % probability of being connected                
                A_prob(graph_element(x,y,m),graph_element(i,j,m))=p; % generate the probability matrix                                        
                        
                if p_random<=p
                    A(graph_element(x,y,m),graph_element(i,j,m))=1;       % connections are symmetric
                    A(graph_element(i,j,m),graph_element(x,y,m))=1;   
                elseif p_random>p
                    A(graph_element(x,y,m),graph_element(i,j,m))=0;       % no connections are symmetric too
                    A(graph_element(i,j,m),graph_element(x,y,m))=0;
                end
            
                %{
                 give the same result :(
                % right border elements            
                p_random=rand;
                p=p0*exp(-(x-(m+1-i))^2/sigma_x -(y-j)^2/sigma_y);           % probability of being connected, other side
                A_prob(graph_element(x,y,m),graph_element((m+1-i),j,m))=p;   % generate the probability matrix                                        
            
                if p_random<=p
                    A(graph_element(x,y,m),graph_element((m+1-i),j,m))=1;       % connections are symmetric
                    A(graph_element((m+1-i),j,m),graph_element(x,y,m))=1;   
                elseif p_random>p
                    A(graph_element(x,y,m),graph_element((m+1-i),j,m))=0;       % no connections are symmetric
                    A(graph_element((m+1-i),j,m),graph_element(x,y,m))=0;
                end                                      
                
               %}                              
            end
        end
    end
                
end


% remove autapses
for i=1:1:m*n
   A(i,i)=0;
%  A_prob(i,i)=0;
end

end

