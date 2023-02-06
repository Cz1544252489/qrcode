tic
%%  part 1
beta_type = 6;

%%   figure part 1
charset={'a','b','c','d'};
filename=cell(1,8);
file={strcat(['p1beta',num2str(beta_type),'fig1','.xls']),strcat(['p1beta',num2str(beta_type),'fig2','.xls'])};
for i=1:8
    filename{i}=strcat(['fig',num2str(floor((i-1)/4)+1),charset{mod((i-1),4)+1},'.xls']);
end

% r=floor(11/30*n_1*n_2*p/(n_1+n_2));

SET=[	1000 1000 36 0.2 1 1 1;	    %(1a)
        1000 1000 36 0.2 0 1 1;	    %(1b)
        1000 1000 36 0.2 1 1 0;	    %(1c)
        1000 1000 36 0.2 0 1 0; 	%(1d)
        1500 1000 44 0.2 1 1 1;     %(2a)
        1500 1000 44 0.2 0 1 1;     %(2b)
        1500 1000 44 0.2 1 1 0;     %(2c)
        1500 1000 44 0.2 0 1 0;     %(2d)
        ];
for i=1:8
    maincode(SET(i,:),file{floor((i-1)/4)+1},charset{mod((i-1),4)+1},beta_type);
    fprintf("p1,%d\n",i);
end

toc

%%   figure part 2
charset={'a','b','c','d'};
filename=cell(1,8);
file={strcat(['p2beta',num2str(beta_type),'fig1','.xls']),strcat(['p2beta',num2str(beta_type),'fig2','.xls'])};
for i=1:8
    filename{i}=strcat(['fig',num2str(floor((i-1)/4)+1),charset{mod((i-1),4)+1},'.xls']);
end

% r=floor(11/30*n_1*n_2*p/(n_1+n_2));

SET=[	2000 2000 73 0.2 1 1 1;	    %(1a)
        2000 2000 73 0.2 0 1 1;	    %(1b)
        2000 2000 73 0.2 1 1 0;	    %(1c)
        2000 2000 73 0.2 0 1 0; 	%(1d)
        2000 1500 62 0.2 1 1 1;     %(2a)
        2000 1500 62 0.2 0 1 1;     %(2b)
        2000 1500 62 0.2 1 1 0;     %(2c)
        2000 1500 62 0.2 0 1 0;     %(2d)
        ];
for i=1:8
    maincode(SET(i,:),file{floor((i-1)/4)+1},charset{mod((i-1),4)+1},beta_type);
    fprintf("p2,%d\n",i);
end

toc

%%   figure part 3
charset={'a','b','c','d'};
filename=cell(1,8);
file={strcat(['p3beta',num2str(beta_type),'fig1','.xls']),strcat(['p3beta',num2str(beta_type),'fig2','.xls'])};
for i=1:8
    filename{i}=strcat(['fig',num2str(floor((i-1)/4)+1),charset{mod((i-1),4)+1},'.xls']);
end

% r=floor(11/30*n_1*n_2*p/(n_1+n_2));
 
SET=[	1000 1000 9 0.05 1 1 1;	    %(1a)
        1000 1000 9 0.05 0 1 1;	    %(1b)
        1000 1000 9 0.05 1 1 0;	    %(1c)
        1000 1000 9 0.05 0 1 0; 	%(1d)
        1500 1000 11 0.05 1 1 1;     %(2a)
        1500 1000 11 0.05 0 1 1;     %(2b)
        1500 1000 11 0.05 1 1 0;     %(2c)
        1500 1000 11 0.05 0 1 0;     %(2d)
        ];
for i=1:8
    maincode(SET(i,:),file{floor((i-1)/4)+1},charset{mod((i-1),4)+1},beta_type);
    fprintf("p3,%d\n",i);
end

toc

%%   figure part 4
charset={'a','b','c','d'};
filename=cell(1,8);
file={strcat(['p4beta',num2str(beta_type),'fig1','.xls']),strcat(['p4beta',num2str(beta_type),'fig2','.xls'])};
for i=1:8
    filename{i}=strcat(['fig',num2str(floor((i-1)/4)+1),charset{mod((i-1),4)+1},'.xls']);
end

% r=floor(11/30*n_1*n_2*p/(n_1+n_2));

SET=[	2000 2000 18 0.05 1 1 1;	    %(1a)
        2000 2000 18 0.05 0 1 1;	    %(1b)
        2000 2000 18 0.05 1 1 0;	    %(1c)
        2000 2000 18 0.05 0 1 0; 	    %(1d)
        2000 1500 15 0.05 1 1 1;        %(2a)
        2000 1500 15 0.05 0 1 1;        %(2b)
        2000 1500 15 0.05 1 1 0;        %(2c)
        2000 1500 15 0.05 0 1 0;        %(2d)
        ];
for i=1:8
    maincode(SET(i,:),file{floor((i-1)/4)+1},charset{mod((i-1),4)+1},beta_type);
    fprintf("p4,%d\n",i);
end

toc

function maincode(SET,filename,sheet,beta_type)

outputflag=0;
if nargin < 2
    filename = 'table.mat';
    if nargin<1
        %     n_1  n_2  r  p   primarydata  initial method way_fsk
        SET=[1000 1000 36 0.2      1           1       1      0   ];
    end
end
%              max_iter epsilon delta  theta
SET_constant = [250     1e-10    1e-4   0.01];
global n_1;
global n_2;
global r;
global data;
global delta;
global Omega;
global theta;
global orth_type;

n_1 = SET(1); n_2 = SET(2); r=SET(3);
p = SET(4);%rates of Omega
M=rankmanQR(n_1,n_2,r);
max_iter=SET_constant(1);%max iteration times
epsilon=SET_constant(2);%stop parameter
delta=SET_constant(3);
theta=SET_constant(4);%parameter to justify weather it owns good orthogonality or not

%%
% orth_type=1 %if we use QR method {1 use 0 don`t}
% way_fsk=0;%apply which line search {1  exact  0 inexact}
% method=0;%apply which descent direction  {1  negetive gradient  0  conjugate direction}
%%
ra=RandStream('mt19937ar','Seed',sum(100*clock));
RandStream.setGlobalStream(ra);
data_=M.gSetData(SET(5));
% generate Omega
Omega=M.gOmega(p);
data=data_.*Omega;%get new data  E(data)=0.01
% fprintf('E(M)=%f\n',sum(sum(data))/nnz(data));
% fprintf('OS=%f\n',nnz(Omega)/M.dim());
settings={  'E(M)',sum(sum(data))/nnz(data);
            'OS',nnz(Omega)/M.dim();
            'n1',n_1;'n2',n_2;
            'r',r;'p',p;
            'max_iter',250;
            'epsilon',1e-10;
            'delta',1e-4;
            'theta',0.01;};
writecell(settings,filename,'Sheet','set','Range','A1');
intial=SET(6);   %1 means spectral initialization  0 means random initial
%beta_type = 1;
if(intial)
    X0=M.gX0(data);
else
    X0=M.rands();
end

method = SET(7);
% way_fsk = SET(8);

for way_fsk = [1,0]
    for orth_type = [1,0]
        %% some important data record
        f0=zeros(1,max_iter);
        s=zeros(1,max_iter);
        RMSE=zeros(1,max_iter);

        X=X0;
        t=1;
        time0=clock;
        if(~method)
            gl=M.grad(X);
            dl=zeros(n_1+n_2,r);
        end
        while(t<=max_iter)
            %function value
            f0(t)=M.f(X);
            %compute graident
            g=M.grad(X);
            if(method)
                d=-g;
            else
                beta=M.fbeta(X,gl,g,dl,beta_type);
                d=-g+beta*dl;
            end
            %find stepsize
            if(way_fsk)
                s(t)=M.fsk(X,-d);
            else
                s(t)=M.ufsk(X,g,d);
            end
            %update iterate
            X=M.next(X,s(t),d);
            %calculate RMSE
            H=Omega.*(X(1:n_1,:)*X(n_1+1:end,:)'-data);
            RMSE(t)=sqrt(sum(sum(H.^2))/nnz(Omega));
            %stop parameter
            %if stop or not
            if(RMSE(t)<epsilon)
                break;
            end
            %update iterate parameter
            t=t+1;
            if(~method)
                dl=d;
                gl=g;
            end
        end
        time1=clock;
        time=(time1-time0)*[0;0;43200;3600;60;1];

        %% output
        temp=RMSE;
        temp(temp==0)=[];
        count=length(temp);
        x=zeros(1,length(RMSE));
        for i=1:count
            x(i)=time/count*i;
        end
        if outputflag==0
            writematrix([x' RMSE'],filename,'Sheet',sheet,'Range','A1');
        end
        if outputflag==1
            writematrix([x' RMSE'],filename,'Sheet',sheet,'Range','C1');
        end
        if outputflag==2
            writematrix([x' RMSE'],filename,'Sheet',sheet,'Range','E1');
        end
        if outputflag==3
            writematrix([x' RMSE'],filename,'Sheet',sheet,'Range','G1');
        end
        outputflag=outputflag+1;

    end
end
end

function [ M ] = rankmanQR( n_1 , n_2 , r )
%   This funtion divides nine parts and the subfuntion in these:
%       1, common part: name, dim, size0
%       2, objective function: f, metric
%       3, metric and gradient: grad
%       4, QR decomposition: qr
%       5, randoms: rands, gX0, gOmega, gSetData
%       6, find stepsize: fsk, ufsk
%       7, update iteration point: next
%       8, conjugate method: fbeta
%       9, output the figure: outplot
    %%   common part
    M.name = @() sprintf('Manifold of %d x %d matrices of rank %d', n_1, n_2, r);
    
    M.dim = @() (n_1+n_2-r)*r;
    
    M.size0 = [n_1,n_2,r];
     
    %% objective function
    M.f=@f;
    function z=f(X)
        global data;
        global Omega;
        s=Omega.*(X(1:n_1,:)*X(n_1+1:n_1+n_2,:)'-data);
        z=0.5*trace(s'*s);
    end
    %% metric and gradient
    %metric
    M.metric = @metric;
    function z =metric( A , B , X )
        global delta;
        global orth_type;
        if nargin <3
            X=zeros(n_1+n_2,r);
        end
        Q=X(1:n_1,:);
        R=X(n_1+1:end,:)';
        if(orth_type)
            z=trace(A(1:n_1,:)'*B(1:n_1,:)*(R*R'+delta*eye(r)))+...
                trace(A(n_1+1:end,:)'*B(n_1+1:end,:)*(1+delta));
        else
            z=trace(A(1:n_1,:)'*B(1:n_1,:)*(R*R'+delta*eye(r)))+...
                trace(A(n_1+1:end,:)'*B(n_1+1:end,:)*(Q'*Q+delta*eye(r)));
        end
    end
    
    %gradient
    M.grad = @grad;
    function Z=grad( X )
        global Omega;
        global data;
        global delta;
        global orth_type;
        Q=X(1:n_1,:);R=X(n_1+1:end,:)';
        S=Omega.*(Q*R-data);
        Z=[S*R';S'*Q];
        if(orth_type)
            Z=[Z(1:n_1,:)/(R*R'+delta*eye(r));Z(n_1+1:end,:)/(1+delta)];
        else
            Z=[Z(1:n_1,:)/(R*R'+delta*eye(r));Z(n_1+1:end,:)/(Q'*Q+delta*eye(r))]; 
        end
    end
    
    %% QR decomposition
    % QR decomposition by MGS
    M.qr = @qr_MGS;
    function [Q,R]=qr_MGS(A)
        [s1,s2]=size(A);
        Q=zeros(s1,r);
        R=zeros(r,s2);
        for k=1:r
            R(k,k)=sqrt(A(:,k)'*A(:,k));
            Q(:,k)=A(:,k)/R(k,k);
            if (k<s2)
                for j=k+1:s2
                    R(k,j)=Q(:,k)'*A(:,j);
                    A(:,j)=A(:,j)-Q(:,k)*R(k,j);
                end
            end
        end
    end
    %% randoms
    % random initialization
    M.rands = @rands;
    function Z=rands(l)
        if nargin<1
            l=0;
        end
        temp=rand(n_1,r+l)*rand(r+l,n_2);
        [Q,R]=M.qr(temp);
        Z=[Q(:,1:r);R(1:r,:)'];
    end

    % generate X for figure process
    M.gXf = @gXf;
    function Z=gXf(data)
        [U,S,V]=svd(data);
        Z=[U(:,1:r)*sqrt(S(1:r,1:r));V(:,1:r)*(sqrt(S(1:r,1:r)))'];
    end

    M.randf = @randf;
    function Z=randf()
        temp=rand(n_1,r)*rand(r,n_2);
        temp=temp/max(temp(:))*255;
        [Q,R]=M.qr(temp);
        Z=[Q(:,1:r);R(1:r,:)'];
    end

    % spectral initialization
    M.gX0 = @gX0;
    function Z=gX0(data)
        [U,S,V]=svd(data);
        Z=[U(:,1:r)*sqrt(S(1:r,1:r));V(:,1:r)*(sqrt(S(1:r,1:r)))'];
    end
    
    % generate Omega
    M.gOmega = @gOmega;
    function Z=gOmega(p)
        Z=rand(n_1,n_2);
        Z(Z>1-p)=1;
        Z(Z<=1-p)=0;
    end
    
    % primary data generation
    M.gSetData = @gSetData;
    function Z=gSetData(input)
        l=0;
        Z=rand(n_1,r+l)*rand(r+l,n_2);  
        if(input==1)
            Z=Z/(max(Z(:))+max(-Z(:)));
        end
        if(input==2)

        end
    end
    
    %% find stepsize
    % by exact line search
    M.fsk=@fsk;
    function [z,ss]=fsk(X,g)
       global Omega;
       global data;
       Q=X(1:n_1,:);R=X(n_1+1:n_1+n_2,:)';
       g_Q=g(1:n_1,:);g_R=g(n_1+1:n_1+n_2,:)';
       A=Omega.*(g_Q*g_R);B=Omega.*(Q*g_R+g_Q*R);O=Omega.*(Q*R-data);
       ss=[2*trace(A'*A) 3*trace(A'*B) 2*trace(A'*O)+trace(B'*B) trace(B'*O)];
       z=min(abs(roots(ss)));
    end
    
    % by inexact line search
    M.ufsk=@ufsk;
    function [z,im]=ufsk(X,g,d)
        %z=1;
        z=M.metric(-g,d,X)/M.metric(d,d,X);
        im=0;im_max=20;
        while(1)
            z=z*0.5^im;
            if(M.f(X)-M.f(X+z*d)+0.001*z*M.metric(g,d,X)>=0)
                break;
            else
                im=im+1;
            end
            if(im>=im_max)
                break;
            end
        end
    end
    
    %% update iteration point
    M.next=@next;
    function [Z]=next(X,s,d)
        global theta;
        global orth_type;
        if isempty(theta)
            disp('don`t set theta default 0.01 '); 
            theta=0.01;
        end
        Z=X+s*d;
        if(orth_type)
            h=abs(trace(Z(1:n_1,:)'*Z(1:n_1,:))-r)/r;
            if(h>theta)
                [Xq,Xr]=M.qr(Z(1:n_1,:));
                Z=[Xq;Z(n_1+1:end,:)*Xr'];
            end
        end
    end

    %% conjugate method
    % find conjugate parameter beta
    M.fbeta=@fbeta;
    function z=fbeta(X,gl,g,dl,type) % gl denotes gradient of last��dl denotes direction of last
        upper_border=1e3;lower_border=1e-3;
        switch type
            case 1  %HS+
                z=min(upper_border,max(lower_border,M.metric(g-gl,g,X)/M.metric(g-gl,dl,X)));
            case 2  %FP
                z=min(upper_border,max(lower_border,M.metric(g,g,X)/M.metric(gl,gl,X)));
            case 3  %PRP
                z=min(upper_border,max(lower_border,M.metric(g-gl,g,X)/M.metric(gl,gl,X)));
            case 4  %CD
                z=min(upper_border,max(lower_border,M.metric(g,g,X)/M.metric(-gl,dl,X)));
            case 5  %LS
                z=min(upper_border,max(lower_border,M.metric(g-gl,g,X)/M.metric(-gl,dl,X)));
            case 6  %DY
                z=min(upper_border,max(lower_border,M.metric(g,g,X)/M.metric(gl,g-gl,X)));
            case 7  %RMIL
                z=min(upper_border,max(lower_border,M.metric(g,g-gl,X)/M.metric(dl,dl,X)));
            case 8  %MRMIL
                z=min(upper_border,max(lower_border,M.metric(g,g-gl-dl,X)/M.metric(dl,dl,X)));
        end
    end
end

