%function that computes the basic measures out of the absorbing Markov
%chain model
%Gregory Roth
%10/10/2016

%Given:
%
%   - U                        a transient matrix
%   - B                        a subset of the transient states
%
%Returns:
%
    %Longevity
  
%   - out.eta                   expected longevity
%   - out.vareta                variance of longevity
%   - out.nuMatrix              mean of the time spent in each state
%   - out.varMatrix             variance of time spent in each state

    %Occupancy time in B
    
%   - out.tau                 mean 
%   - out.tau2                second moment 
%   - out.vartau              variance 
%   - out.cv_tau                  coefficient of variation 

    %Reaching the subset B

%   - out.p_a                   probability to reach
%   - out.ptilde_a              extension of p_a
%   - out.t_B                   mean time to reach 
%   - out.vart_B                variance time ot reach

    %Returning to the subset B

%   - out.p_r                   probability to return
%   - out.lambda                mean time to return
%   - out.varlambda             variance time ot return

    %Useful matrices

%   - out.U_B                   transient matrix of the sub Markov cahin
%   - out.U_c
%   - out.W_B
%   - out.W_Bc
%   - out.Q
%   - out.L
%   - out.K;
%   - out.Atilde_k
%   - out.q                     position-in-the-rearranged-matrix vector

%***********************************************************************

function out=abs_MC(U,B)
siz=size(U);
s=siz(1);                   %size of U

sa=length(B);               %number of states in the subset B
st=s-sa;                    %number of states in the complement of B


%useful vectors
es=ones(1,s);
esa=ones(1,sa);
est=ones(1,st);

%permutation vector, used to rearrange the matrix U
p=(1:s);
for i=1:sa
    p(p==B(i))=[];
end
p=[p,B];

%Rearranging the indices of the transient states such that the last states are the states in B
        Utemp=U;
        for i=1:s
            Utemp(:,i)=U(:,p(i));
        end
        Utemp2=Utemp;
        for i=1:s
        Utemp2(i,:)=Utemp(p(i),:);
        end
        Uprime=Utemp2; %rearranged transition matrix

        
       
%Decomposing the transition matrix and computing useful matrices

        U_k=Uprime(1:st,1:st);  %transitions from B^c to B^c
        K=Uprime(st+1:end,1:st);%transitions from B^c to B
        L=Uprime(1:st,st+1:end);%transitions from B to B^c
        Q=Uprime(st+1:end,st+1:end);%transitions from B to B
        
        U_sub=Uprime(:,st+1:end); %transition probabilities from states in B to any states (except death)
        
        Nprime=inv(speye(s)-Uprime); %fundamental matrix of the original chain (rearranged)
        N=inv(speye(s)-U); %fundamental matrix of the original chain (non rearranged)
        
        N_k=inv(eye(st)-U_k); %fundamental matrix for of the killed MC

%Absorbtion probabilities, via the killed MC

        A_k=K*N_k;
        Atilde_k=[A_k,eye(sa)]; %extension of the absorbing probability matrix
        
         p_a=esa*A_k;     %probabilities of abosrbition by B = reaching B probabilities
         ptilde_a=esa*Atilde_k; %extension of p_a

%creating the conditonal Markov chain  
         D_a=diag(p_a);
            
        %transient matrix and mortality matrix
         invpa=p_a.^(-1);
         invpa(isinf(invpa)) = 0; %replace NaN by zero
         invDa=diag(invpa);
         U_c=D_a*U_k*invDa;
         c=(esa*K)/D_a;
         
         %fundamental matrix
         N_c=inv(speye(st)-U_c);

%creating the sub Markov chain 

        %transient matrix
        U_B=Atilde_k*U_sub;

        %fundamental matrix
        N_B=inv(speye(sa)-U_B);
        
%computing the measures

        %longevity
        nuMatrix=N;
        eta=sum(N);
        vareta=es*N*(2*N-eye(s))-eta.*eta;
        varMatrix=(2*diag(diag(N))-eye(s))*N-N.*N;
        
        %Occupancy time in B
            
        tau_B=sum(N_B);
        tau2_B=esa*N_B*(2*N_B-eye(sa));
        nuMatrix_Q=N_B;
        varMatrix_Q=(2*diag(diag(N_B))-eye(sa))*N_B-N_B.*N_B;
        
        tau=tau_B*Atilde_k;
        tau2=tau2_B*Atilde_k;
        vartau=tau2-tau.*tau;

        %Reaching the subset B

        t_B=est*N_c;        %mean time to reach
        t2_B=est*N_c*(2*N_c-eye(st));
        vart_B=t2_B-t_B.*t_B;   

        %Returning to the subset B

        p_r=esa*U_B;  %return probabilities 
        D_r=diag(p_r);
        
        W_Bc=(D_a*L)/D_r;   %conditional transition probabilities B to B^c given individual returns in B
        W_B=(Q)/D_r;   %conditional transition probabilities B to B given individual returns in B
        
        lambda=esa+t_B*W_Bc;    %mean time to return
        varlambda=t2_B*W_Bc-(t_B*W_Bc).*(t_B*W_Bc); %variance time ot return


%rearrange the output vectors in the initial order

        q=(1:s);
        for i=1:s
            q(i)=find(p==i);
        end
        
        tautemp=tau;
        tau2temp=tau2;
        vartautemp=vartau;

        for i=1:s
            tautemp(i)=tau(q(i));
            tau2temp(i)=tau2(q(i));
            vartautemp(i)=vartau(q(i));
        end

%outputs

out.eta=eta;                   
out.vareta =vareta;             
out.nuMatrix= nuMatrix;              
out.varMatrix= varMatrix;             

out.tau=tautemp;                
out.tau2=tau2temp;                
out.vartau=vartautemp;              
out.cv_tau= sqrt(vartautemp)./tautemp;                  

out.p_a=p_a;
out.ptilde_a= ptilde_a;
out.t_B=t_B;                   
out.vart_B=vart_B;               

out.p_r=p_r;                   
out.lambda=lambda;              
out.varlambda=varlambda;   

out.U_B=U_B;
out.U_c=U_c;
out.W_B=W_B;
out.W_Bc=W_Bc;
out.Q=Q;
out.L=L;
out.K=K;
out.Atilde_k=Atilde_k;
out.q=q;
end


