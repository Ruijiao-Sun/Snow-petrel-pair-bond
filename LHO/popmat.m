%absorbing markov chain to calculate LHO of snow petrel
%theta: parameter estimates
%current population structure
function [U,I,P,N,R1,R2,R3,R4] = popmat(theta)
%%survival
%S         NW        ND     NB      W          D
S=[theta(1)    0         0       0      0          0                 
     0      theta(2)     0       0      0          0                   
     0         0     theta(3)    0      0          0                 
     0         0         0     theta(4) 0          0                  
     0         0         0       0    theta(5)     0                  
     0         0         0       0      0     theta(6)];

%%widowhood
W=[ 1-theta(7)    0          0       0         0     0
          0    1-theta(8)    0       0         0     0     
          0       0      1-theta(9)  0         0     0         
          0       0          0     1-theta(10) 0     0      
          0       0          0       0         1     0      
          0       0          0       0         0     1  
     theta(7)     0          0       0         0     0     
          0     theta(8)     0       0         0     0    
          0       0      theta(9)    0         0     0    
          0       0          0    theta(10)    0     0    ];

%%divorce
D=[1-theta(13)   0         0         0         0    0    0    0    0    0
      0     1-theta(14)    0         0         0    0    0    0    0    0
      0          0    1-theta(15)    0         0    0    0    0    0    0
      0          0         0    1-theta(16)    0    0    0    0    0    0
      0          0         0         0         1    0    0    0    0    0
      0          0         0         0         0    1    0    0    0    0
      0          0         0         0         0    0    1    0    0    0
      0          0         0         0         0    0    0    1    0    0
      0          0         0         0         0    0    0    0    1    0
      0          0         0         0         0    0    0    0    0    1
     theta(13)   0         0         0         0    0    0    0    0    0
      0       theta(14)    0         0         0    0    0    0    0    0
      0          0     theta(15)     0         0    0    0    0    0    0
      0          0         0     theta(16)     0    0    0    0    0    0];


%%Breeding
B=[theta(19)  theta(20)  theta(21)  theta(22)    0         0          0         0         0          0        0         0         0          0   
      0         0           0          0       theta(23)   0      theta(19) theta(20) theta(21)  theta(22)    0         0         0          0
      0         0           0          0         0      theta(24)     0         0         0          0      theta(19) theta(20) theta(21)  theta(22)
 1-theta(19) 1-theta(20) 1-theta(21) 1-theta(22) 0         0          0         0         0          0        0         0         0          0
      0         0           0          0      1-theta(23)  0  1-theta(19) 1-theta(20) 1-theta(21) 1-theta(22) 0         0         0          0
      0         0           0          0         0     1-theta(24)    0         0         0          0  1-theta(19) 1-theta(20) 1-theta(21) 1-theta(22) ];
 

U = (S'*W'*D'*B')';

f1 = theta(31:33)*U(1:3,1); %number of offspring,  not real fertility, real fertility need to time 0.5
f2 = theta(31:33)*U(1:3,2);
f3 = theta(31:33)*U(1:3,3);
f4 = theta(31:33)*U(1:3,4);
f5 = theta(31:33)*U(1:3,5);
f6 = theta(31:33)*U(1:3,6);

%A = U + F; 
 m = 1-sum(U);
  [s,s] = size(U); 
      I = eye(s);
      P = [U zeros(s,1);m 1];      %absorbing markov chain
      N = inv(I-U); %fundamental matrices
      
%r1 = theta(31:36);R1=[ones(s+1,1).*r1 zeros(s+1,1)]; %depends on
r1 = [f1 f2 f3 f4 f5 f6]; %vector of fertility
                         R1=[ones(s+1,1).*r1 zeros(s+1,1)];
r2bern=r1;               R2=[ones(s+1,1).*r2bern zeros(s+1,1)];
r3bern=r1;               R3=[ones(s+1,1).*r3bern zeros(s+1,1)];
r4bern=r1;               R4=[ones(s+1,1).*r4bern zeros(s+1,1)];
end

