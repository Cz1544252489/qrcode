ccc;% clc,clear,close
%%
% This code is builded in 2022/12/1 for ManQR to promote performance
% And last debugging or promoting time is 2022/12/6
% By Chenzhuo

%% alternative 1 for Synthetic data
% S for stable setting on iteration, D for data and F for function sets
%
n1 = 1000;    % row number
n2 = 1000;    % column number
r = 36;       % rank of Synthetic matrix
p = 0.2;      % frequency of Synthetic matrix
max_iter = 250;   % max iteration times
epsilon = 1e-10;  % stopping criterion
delta = 1e-4;     % metric needs
theta = 0.01;     % if need to process by QR
S = struct("n1",n1,"n2",n2,"r",r,"p",p,"max_iter",max_iter,...
    "epsilon",epsilon,"delta",delta,"theta",theta);

% build 'observed' data
RaStr=RandStream('mt19937ar','Seed',sum(100*clock)); % according to time
RandStream.setGlobalStream(RaStr);         % choose random starting point
M = rand(S.n1,S.r)*rand(S.r,S.n2);      % build Synthetic matrix
if (true), M=M/max(M,[],'all');  end    % change the range of M if needs
D.Omega = rand(S.n1,S.n2) < S.p;        % build Omega matrix for observing
D.PM = M.*(D.Omega);                    % Input matrix P_Omega(M)

%% alternative 2 for real data or image


%% before iteration
% C for changeable setting on iteration
% C.ot = 1;        % 1 for using orthogonality and 0 for not
% C.dm = "cd";     % {"sd","cd"} sd for steepest dircetion and cd for conjugate direction
C.dmb = "LS";    % beta choice {"HS+","FP","PRP","CD","LS","DY","RMIL","MRMIL"} in conjugate method
% C.ls = "iels";    % {"els","iels"}
C.in = "sv";     % {"sv","ra"}  svdPM or random intialization

setting = ["sd","els";
    "sd","iels";
    "cd","els";
    "cd","iels";];
%% iteration part
% color of line and point
P=struct('line','r','point','b');
for j = 1:4
    C.dm = setting(j,1);
    C.ls = setting(j,2);
    for i=[1 0]

        C.ot = i;

        F = rankmanQR(C,S,D);     % load function sets F

        [RMSE, time] = iteration(C, S, D, F);

        % plotfig(filename,P);

    end
end

%% function sets
function [RMSE, time] = iteration(C, S, D, F)

% initialization
switch C.in
    case "ra"
        X = F.ra();
    case "sv"
        X = F.sv();
end

% for same start point on figure
RMSE=zeros(1,S.max_iter);
H = (D.Omega).*(X(1:S.n1,:)*X(S.n1+1:S.n1+S.n2,:)'-D.PM);
RMSE(1) = sqrt(sum(H.^2,"all")/nnz(D.Omega));

time0 = clock;

if(C.dm=="cd")
    gl = F.g(X);
    dl = zeros(S.n1+S.n2,S.r);
end

% start from 2
for t=2:S.max_iter

    % function value
    % f0=F.f(X);

    % compute gradient
    g = F.g(X);

    % descent method
    switch C.dm
        case "sd"
            d = -g; % steepest
        case "cd"
            beta = F.fb(X,gl,g,dl);
            d = -g+beta*dl; % conjugate
    end

    % line search
    switch C.ls
        case "els"  % exact line search
            s = F.els(X,-d);
        case "iels" % inexact line search
            s = F.iels(X,g,d);
    end

    % update iterate
    X = F.n(X,s,d);

    % calculate RMSE
    H = (D.Omega).*(X(1:S.n1,:)*X(S.n1+1:S.n1+S.n2,:)'-D.PM);
    RMSE(t) = sqrt(sum(H.^2,"all")/nnz(D.Omega));
    disp(RMSE(t));
    if RMSE(t) < S.epsilon, break; end

    % update iteration parameter
    if(C.dm=="cd")
        dl = d;
        gl = g;
    end
end

rebuildedM = X(1:S.n1,:)*X(S.n1+1:S.n1+S.n2,:)';

time1 = clock;
time=(time1-time0)*[0;0;43200;3600;60;1];

% save data settings and result
dt = char(datetime);
char0 = regexp(dt,'(-|\s|:)');
for i=1:length(char0),dt(char0(i))='_';end
dt = string(dt);
saveplace='L:\Algorithm_data\';
filename = dt+'_or'+string(C.ot)+'_'+C.dm+'_'+C.dmb+'_'+C.ls+'_'+C.in+'_'+...
    string(S.n1)+'_'+string(S.n2)+'_r'+string(S.r)+'_p'+string(100*S.p)+'.mat';
save(saveplace+filename,"RMSE","time","dt","D","S","C","rebuildedM",'-mat');

end

function F = rankmanQR(C,S,D)

% random initialization
F.ra = @rands;
    function z = rands()
        Tmp = rand(S.n1,S.r)*rand(S.r,S.n2);
        [Q,R] = F.qr(Tmp);
        z = [Q(:,1:S.r);R(1:S.r,:)'];
    end

% svd PM initialization
F.sv = @svdPMInitiailization;
    function z = svdPMInitiailization()
        [U,Si,V]=svd(D.PM);
        z = [U(:,1:S.r)*sqrt(Si(1:S.r,1:S.r));V(:,1:S.r)*(sqrt(Si(1:S.r,1:S.r)))'];
    end

% objective function
F.f=@f;
    function z = f(X)
        s = (D.Omega).*(X(1:S.n1,:)*X(S.n1+1:S.n1+S.n2,:)'-D.PM);
        z = 0.5*norm(s,'fro')^2;
    end

% metric
F.m = @metric;
    function z = metric(A,B,X)
        if nargin<3
            X = zeros(S.n1+S.n2,S.r);
        end
        Q = X(1:S.n1,:);R = X(S.n1+1:S.n1+S.n2,:)';
        if(C.ot)
            z=trace(A(1:S.n1,:)'*B(1:S.n1,:)*(R*R'+(S.delta)*eye(S.r)))+...
                trace(A(S.n1+1:end,:)'*B(S.n1+1:end,:)*(1+S.delta));
        else
            z=trace(A(1:S.n1,:)'*B(1:S.n1,:)*(R*R'+(S.delta)*eye(S.r)))+...
                trace(A(S.n1+1:end,:)'*B(S.n1+1:end,:)*(Q'*Q+(S.delta)*eye(S.r)));
        end
    end

% gradient
F.g = @gradient;
    function z = gradient(X)
        Q = X(1:S.n1,:);R = X(S.n1+1:S.n1+S.n2,:)';
        Tmp = (D.Omega).*(Q*R-(D.PM));
        z=[Tmp*R';Tmp'*Q];
        if(C.ot)
            z=[z(1:S.n1,:)/(R*R'+(S.delta)*eye(S.r));z(S.n1+1:end,:)/(1+S.delta)];
        else
            z=[z(1:S.n1,:)/(R*R'+(S.delta)*eye(S.r));z(S.n1+1:end,:)/(Q'*Q+(S.delta)*eye(S.r))];
        end
    end

% QR decomposition by MGS
F.qr = @qr_MGS;
    function [Q,R]=qr_MGS(A)
        [s1,s2]=size(A);
        Q=zeros(s1,S.r);
        R=zeros(S.r,s2);
        for k=1:S.r
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

% find stepsize by exact line search
F.els = @exact_line_search;
    function z = exact_line_search(X,g)
        Q = X(1:S.n1,:);
        R = X(S.n1+1:S.n1+S.n2,:)';
        gQ = g(1:S.n1,:);
        gR = g(S.n1+1:S.n1+S.n2,:)';
        A = (D.Omega).*(gQ*gR);
        B = (D.Omega).*(Q*gR+gQ*R);
        O = (D.Omega).*(Q*R-D.PM);
        ss = [2*norm(A,"fro")^2 3*trace(A'*B) 2*trace(A'*O)+norm(B,"fro")^2 trace(B'*O)];
        z=min(abs(roots(ss)));
    end

% find stepsize by exact line search
F.iels = @inexact_line_search;
    function z =inexact_line_search(X,g,d)
        z = F.m(-g,d,X)/F.m(d,d,X);
        im_max = 20;
        for im=0:im_max
            z = z*0.5^im;
            if(F.f(X)-F.f(X+z*d)+0.001*z*F.m(g,d,X)>=0)
                break;
            end
        end
    end

% update iteration point
F.n = @next_iteration_point;
    function z = next_iteration_point(X,s,d)
        z = X+s*d;
        if(C.ot)
            h = abs(norm(z(1:S.n1,:),"fro")^2/S.r-1);
            if(h>S.theta)
                [Xq,Xr] = F.qr(z(1:S.n1,:));
                z = [Xq;z(S.n1+1:S.n1+S.n2,:)*Xr'];
            end
        end
    end

% conjugate method
F.fb = @find_beta_in_conjugate_method;
    function z = find_beta_in_conjugate_method(X,gl,g,dl)
        % gl denotes last gradient and dl denotes last direction
        max_value = 1e3;
        min_value = 1e-3;
        switch C.dmb
            case "HS+"
                z = F.m(g-gl,g,X)/F.m(g-gl,dl,X);
            case "FP"
                z = F.m(g,g,X)/F.m(gl,gl,X);
            case "PRP"
                z = F.m(g-gl,g,X)/F.m(gl,gl,X);
            case "CD"
                z = F.m(g,g,X)/F.m(-gl,dl,X);
            case "LS"
                z = F.m(g-gl,g,X)/F.m(-gl,dl,X);
            case "DY"
                z = F.m(g,g,X)/F.m(gl,g-gl,X);
            case "RMIL"
                z = F.m(g,g-gl,X)/F.m(dl,dl,X);
            case "MRMIL"
                z = F.m(g,g-gl-dl,X)/F.m(dl,dl,X);
        end
        z = max(min_value,min(max_value,z));
    end

end

function PM = ImgPro(filename, p)
Img = imread(filename);
subplot(2,4,1);
imshow(Img);
title('real image');
Imggray = rgb2gray(Img);
subplot(2,4,2);
imshow(Imggray);
title('gray image');
Omega = rand(size(Imggray)) < p;
Omega = uint8(Omega);
PM = Imggray.*Omega;
subplot(2,4,3);
imshow(PM);
title("'obseved' image");
PM = double(PM);
end


