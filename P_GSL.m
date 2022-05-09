function mu=P_GSL(Y, fs)
%% The proposed P-GSL method

%% Estimate the period
N_orginal=size(Y,1);
P=1/search_P(Y,fs);
R=fs*P;

L= floor(N_orginal/R);
Ini_pints=[0:L-1]*R+1;
Ini_pints=round(Ini_pints);
R=round(R);
y_new=[];
index_record=[];
for ii=1:length(Ini_pints)
    y_new=[y_new;Y(Ini_pints(ii):Ini_pints(ii)+ R-1 )];
    index_record=[index_record; [Ini_pints(ii):Ini_pints(ii)+ R-1]'  ];
end

%% Normalized
[N,T]=size(y_new);
Norm_y=norm(y_new,'fro')/sqrt(N*T);
y_new=y_new/Norm_y;


%%  Initialization
a=1e-10;
b=a;
gamma_with_tau=ones(N,1);
tau=ones(L,1);
Z=ones(N,3)/3;
alpha=1;
iter = 0;
maxiter=100;
rho=0.5;
converged=false;

while ~converged
    %%  Update x
    deltal=[gamma_with_tau(2:length(gamma_with_tau));gamma_with_tau(1)];
    deltar=[gamma_with_tau(end);gamma_with_tau(1:length(gamma_with_tau)-1)];
    gamma_with_tau=(  Z(:,1).*deltar + Z(:,2).*gamma_with_tau +  Z(:,3).*deltal   );
    Sigma=1./(alpha + gamma_with_tau );
    mu= alpha*( Sigma.* y_new);
    
    %% Update alpha
    resid=y_new-mu;
    alpha=( T*N + 2*a )/( 2*b +  norm(resid, 'fro')^2  + T*  sum(Sigma)  );
    
    %% Update gamma
    sum_mu=sum( mu.*conj(mu), 2);
    mu2=sum_mu + T*real(Sigma);
    mu2l=[mu2(2:length(mu2));mu2(1)];
    mu2r=[mu2(end);mu2(1:length(mu2)-1)];
    Z_2=Z(:,2);
    Z_1l=[Z(2:length(mu2),1);Z(1,1)];
    Z_3r=[Z(end,3);Z(1:length(mu2)-1,3)];
    mu2_combine= Z_1l.*mu2l + Z_2.*mu2 + Z_3r.*mu2r;
    
    mu2_reshape=reshape(mu2_combine, R, L );
    Z_shape=reshape((Z_1l + Z_2 + Z_3r), R, L );
    c_k=( sum(Z_shape,2)*T+  2*a   );
    d_k=(  sum(mu2_reshape*diag(tau),2)  +  2*b  );
    gamma_core=c_k./d_k;
    gamma_with_tau=kron(tau,gamma_core);
    ln_gamma_core= psi( c_k ) -  log(  d_k  );
    ln_delta_with_ones=kron(ones(L,1),ln_gamma_core) ;
    
    %% Update Z
    Z_old=Z;
    ln_deltal=[ln_delta_with_ones(2:length(ln_delta_with_ones));ln_delta_with_ones(1)];
    ln_deltar=[ln_delta_with_ones(end);ln_delta_with_ones(1:length(ln_delta_with_ones)-1)];
    deltal=[gamma_with_tau(2:length(gamma_with_tau));gamma_with_tau(1)];
    deltar=[gamma_with_tau(end);gamma_with_tau(1:length(gamma_with_tau)-1)];
    t1 = ln_deltar - mu2.*deltar ;
    t2 = ln_delta_with_ones  - mu2.*gamma_with_tau ;
    t3 = ln_deltal - mu2.*deltal ;
    et=[t1,t2,t3]/2;
    temp_p= exp(et);
    Z= ((  1./sum(temp_p,2) )*ones(1,3)).*   temp_p;
    Z=rho*Z + (1-rho)*Z_old;  % damping
    
    if iter>=30
        [~,ind_sort]=sort(gamma_core,'descend');
        gamma_core(ind_sort(1:round(end*0.8)))=1e3;
        if iter>=maxiter-5
            gamma_core(ind_sort(1:round(end*0.8)))=1e5;
        end
    end
    
%     figure(10);plot(1./gamma_core)
    
    %% Update tau
    delta_11=kron(ones(L,1),gamma_core);
    deltal=[delta_11(2:length(delta_11));delta_11(1)];
    deltar=[delta_11(end);delta_11(1:length(delta_11)-1)];
    delta_temp=(  Z(:,1).*deltar + Z(:,2).*delta_11 +  Z(:,3).*deltal   );
    temp2=reshape(delta_temp.*mu2, R, L );
    cl= sum(temp2,1)';
    tau_old=tau;
    tau= ( R*T )./( cl );
    tau=rho*tau+ (1-rho)*tau_old; % damping
    gamma_with_tau=kron(tau,gamma_core);
    
    if iter >= maxiter
        converged = true;
    end
    iter = iter + 1;
    
end

mu=mu*Norm_y;
%% calibrate the position
y_res=zeros(N_orginal,T);
for ii=1:length(mu)
    y_res(index_record(ii))=mu(ii) ;
end
mu=y_res;