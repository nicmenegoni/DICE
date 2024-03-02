function [mu_set,mu,I,rho,P21,MTL2D,stdevL2D,max_L_uncut,min_L_uncut,max_l_cut,min_l_cut] = CWstatistic(l_old,L_old,Set,m_mauldon_old,n_mauldon_old,CWr,countCWDisc_inter)
m_mauldon=zeros(length(Set),1);
n_mauldon=zeros(length(Set),1);
L=zeros(length(Set),1);
l=zeros(length(Set),1);
for i =1 : length(Set)
    L(i) = max(L_old(i,:));
    l(i) = max(l_old(i,:));
    for j= 1 :length(L_old(i,:))
        if L_old(i,j)==max(L_old(i,:))
            m_mauldon(i) = m_mauldon_old(i,j);
            n_mauldon(i) = n_mauldon_old(i,j);
        end
    end
end
%% CALULATION PARAMTERS (MTL, St.DEV, P21, mu, I, etc.)
nplane = length(Set);
L_corr=zeros(nplane, max(Set));
l_corr=zeros(nplane, max(Set));

for i=1:nplane
    for j=1:max(Set)
        if Set(i)==j
            if l(i)>0
                l_corr(i,j)=l(i)';
                L_corr(i,j)=L(i)';
            else
                l_corr(i,j)=0;
                L_corr(i,j)=0;
            end
        end
    end
end
m_mauldon2=zeros(nplane, max(Set));
n_mauldon2=zeros(nplane, max(Set));

for i=1:nplane
    for j=1:max(Set)
        if Set(i)==j
            m_mauldon2(i,j)=m_mauldon(i);
            n_mauldon2(i,j)=n_mauldon(i);
        end
    end
end

for j=1:(max(Set))
    M_mauldon(1,j)=sum(m_mauldon2(:,j));
    N_mauldon(1,j)=sum(n_mauldon2(:,j));
end



for j=1:max(Set)
    if (M_mauldon(1,j)+N_mauldon(1,j))>(sum(countCWDisc_inter(:,j)));;
        N_mauldon(1,j)=N_mauldon(1,j)-(M_mauldon(1,j)+N_mauldon(1,j)- sum(countCWDisc_inter(:,j)));
    end
end

% Calculate Mauldon estimator and P21 for the single sets
for j=1:max(Set)
    % mu = Mean Trace Length
    mu_set(1,j)=((pi*CWr)/2)*(N_mauldon(1,j)/M_mauldon(1,j));
    
    % I = trace intesity
    I_set(1,j)=N_mauldon(1,j)/(4*CWr);
    
    % rho = trace density
    rho(1,j)=M_mauldon(1,j)/((pi*CWr)/2);
    P21_set(1,j)=sum(l_corr(:,j))/(pi*CWr^2);
end

mu_set(isnan(mu_set))=0;
% Calculate Mauldon estimator for all the sets together
mu = ((pi*CWr)/2)*(sum(N_mauldon)/sum(M_mauldon));

mu(isnan(mu))=0;

I = sum(N_mauldon)/(4*CWr);
rho = sum(M_mauldon)/((pi*CWr)/2);
P21=sum(sum(l_corr(:,:)))/(pi*CWr^2);



for j=1:max(Set)
    clear Ltemp
    Ltemp=L_corr;
    MTL2D_set(1,j)=mean(Ltemp(Ltemp>0));
    stdev2D_set(1,j)=std(Ltemp(Ltemp>0));
    max_L_uncut_set(1,j)=max(Ltemp(Ltemp>0));
    min_L_uncut_set(1,j)=min(Ltemp(Ltemp>0));
    max_l_cut_set(1,j)=max(Ltemp(Ltemp>0));
    min_l_cut_set(1,j)=min(Ltemp(Ltemp>0));
end
MTL2D=mean(L_corr(L_corr>0));
stdevL2D=std(L_corr(L_corr>0));
max_L_uncut = max(L_corr(L_corr>0));
min_L_uncut = min(L_corr(L_corr>0));
max_l_cut = max(l_corr(l_corr>0));
min_l_cut = min(l_corr(l_corr>0));
end
