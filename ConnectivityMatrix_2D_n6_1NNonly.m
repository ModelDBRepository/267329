%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up connectivity matrix W

% 2D network
% set network size n 
%n = 8;

n=6;

% different weight especially for weak coupling
eps = 1; %this is "d" in our paper
eps1 = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% one neighbor coupling
k = 1;  % set which neighbor
A = zeros(n);
    A = A +diag(ones(n-k,1),k)+diag(ones(n-k,1),-k); %for 2nd NN h/v

% add synapses to end cells = periodic boundary conditions
for m=k:-1:1
        A(n-(k-m),m) = 1;
        A(m,n-(k-m)) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B=A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% connectivity matrix using Kronecker Tensor Product - 1st NN in all directions:
% horizontal & vertical directions + and diagonal direction

% Use W1 for Fig. 5 
% Use W for Fig. 6
W=kron(A,eye(n))+kron(eye(n),B) + eps*kron(A,B); % diagonal coupling
W1=kron(A,eye(n))+kron(eye(n),B) + eps1*kron(A,B); % without diagonal coupling
