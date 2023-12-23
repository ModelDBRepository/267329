%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up connectivity matrix W

% 2D network
% set network size n 
n = 6;

% different weight especially for weak coupling
eps = 1; %turn on parameter for 2nd NN coupling
eps1 = 1; %turn on parameter for 1st NN coupling
eps2= 4;
eps3=0.4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% one neighbor coupling
k = 2;  % set which neighbor
k1 = 1;
A = zeros(n);
A1 = zeros (n);
A2 = zeros (n);
    A = A +diag(ones(n-k,1),k)+diag(ones(n-k,1),-k); %for 2nd NN h/v
    A1 = A1 +diag(ones(n-k1,1),k1)+diag(ones(n-k1,1),-k1); % for 1st NN all directions
    A2 = A2 +diag(ones(n-k1,1),k1)+diag(ones(n-k1,1),-k1); 
    % for 1st NN all directions with heterogeneous coupling strengths

% add synapses to end cells = periodic boundary conditions
for m=k:-1:1
        A(n-(k-m),m) = 1;
        A(m,n-(k-m)) = 1;
end

for m=k1:-1:1
        A1(n-(k1-m),m) = 1;
        A1(m,n-(k1-m)) = 1;
end

for m=k1:-1:1
        A2(n-(k1-m),m) = eps3;
        A2(m,n-(k1-m)) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=A;
B1=A1;
B2=A2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% connectivity matrix using Kronecker Tensor Product - horizontal &
% vertical directions + and diagonal direction
% Use W2 for Fig. 7, Fig. 8A
% Use W4 for Fig. 8B
% Use W5 for Fig. 8C
% with d1>0
W2=eps*(kron(A,eye(n))+kron(eye(n),B)) + eps1*(kron(A1,eye(n))+kron(eye(n),B1)+kron(A1,B1));

% with d1=0
W3=eps*(kron(A,eye(n))+kron(eye(n),B)) + eps1*(kron(A1,eye(n))+kron(eye(n),B1));

% with d1=4=h2=v2 
W4=eps2*(kron(A,eye(n))+kron(eye(n),B)+kron(A1,B1)) + eps1*(kron(A1,eye(n))+kron(eye(n),B1));

% with v1=0.4*h1 and h1=d=h2=v2
W5=eps*(kron(A,eye(n))+kron(eye(n),B)) + eps1*(eps3*kron(A1,eye(n))+kron(eye(n),B1)+kron(A1,B1));
