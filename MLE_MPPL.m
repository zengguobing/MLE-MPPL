function [ phase_series_MLE_MPPL,pol_detR ] = MLE_MPPL( vv_cell,vh_cell,homo_num, homo_num_thresh, homo_index, homo_wnd_rg, homo_wnd_az, n_images )
warning('off') ;
n_pol = 2;%%number of polarimetric channels
nr=size(vh_cell{1,1},1);
nc=size(vh_cell{1,1},2);
u=zeros(n_images*n_pol,1);
u_pol=zeros(n_pol,1);
u_intf=zeros(n_images,1);
phase_series_MLE_MPPL=cell(1,n_images);%% MLE_MPPL phase series to be estimated
%%phase series initialization
for ii=1:n_images
    phase_series_MLE_MPPL{1,ii} = vv_cell{1,ii};
end
pol_detR=zeros(nr,nc);%%pol-detR value
radius_rg=(homo_wnd_rg - 1) / 2;
radius_az=(homo_wnd_az - 1) / 2;

fbar = waitbar(0,'Please wait...');
for ii=1:nr
    for jj=1:nc
        if homo_num(ii, jj) < homo_num_thresh %%processing DS pixels
            continue;
        end
        mask=ones(homo_wnd_az, homo_wnd_rg);
        mask(:,:)=homo_index(ii,jj,:,:);
        Covariance=0;
        Cov_vh=0;Cov_vv=0;
        Cov_pol=0;
        count=0;
        for iii=ii-radius_az:ii+radius_az
            for jjj=jj-radius_rg:jj+radius_rg
                if iii < 1 || iii > nr || jjj < 1 || jjj > nc
                    continue;
                end
                ixx=iii-(ii-radius_az)+1;iyy=jjj-(jj-radius_rg)+1;
                if mask(ixx, iyy)<1
                    continue;
                end
                count=count+1;
                for kk=1:n_images
                    u(kk)=vv_cell{1,kk}(iii,jjj)/sqrt(2);
                end
                for kk=n_images+1:n_images*2
                    u(kk)=vh_cell{1,kk-n_images}(iii,jjj)*sqrt(2);
                end
                Covariance=Covariance+u*u';
                for kk=1:n_images
                    u_intf(kk)=vv_cell{1,kk}(iii,jjj)/sqrt(2);
                end
                Cov_vv=Cov_vv+u_intf*u_intf';
                for kk=1:n_images
                    u_intf(kk)=sqrt(2)*vh_cell{1,kk}(iii,jjj);
                end
                Cov_vh=Cov_vh+u_intf*u_intf';
                
                for kk=1:n_images
                    u_pol(1)=vv_cell{1,kk}(iii,jjj)/sqrt(2);
                    u_pol(2)=vh_cell{1,kk}(iii,jjj)*sqrt(2);
                    Cov_pol=Cov_pol+u_pol*u_pol';            
                end 
            end
        end
        Cov_pol=Cov_pol/n_images/count;
        Covariance=Covariance/count;
        Cov_vh=Cov_vh/count;
        Cov_vv=Cov_vv/count;
        Covariance = sqrt(diag(1./diag(Covariance)))*Covariance*sqrt(diag(1./diag(Covariance)));
        Cov_pol = sqrt(diag(1./diag(Cov_pol)))*Cov_pol*sqrt(diag(1./diag(Cov_pol)));
        Cov_vh = sqrt(diag(1./diag(Cov_vh)))*Cov_vh*sqrt(diag(1./diag(Cov_vh)));
        Cov_vv = sqrt(diag(1./diag(Cov_vv)))*Cov_vv*sqrt(diag(1./diag(Cov_vv)));
        C_coh=0.5*Cov_vv+0.5*Cov_vh;

        %%%MLE-MPPL processing
        [B, C] = SKPDecomposition( Covariance,[n_pol n_pol],[n_images n_images], 1 );
        XX=0;
        XXX=0;
        for mm=1:n_pol * n_pol
            XX=XX+trace(inv(Cov_pol)*B{1,mm})*(inv(abs(C_coh)).*C{1,mm});
            XXX=XXX+trace(inv(Cov_pol)*B{1,mm})*C{1,mm}/n_pol;
        end
        [X,CC]=eig(XX);
        lambda = diag((real(CC)));
        lambda(lambda<0)=1e10;
        [~,lambda_ii] = min(lambda);
        x = (X(:,lambda_ii));
        for kk=1:n_images
            phase_series_MLE_MPPL{1,kk}(ii,jj)=x(kk);
        end
        %%%pol-detR computation
        L1=0.0;
        gamma_MLE=real(diag(conj(exp(1j*angle(x))))*XXX*(diag(conj(exp(1j*angle(x)))))');
        for mm=1:n_pol*n_pol
            temp1=diag(conj(exp(1j*angle(x))))*C{1,mm}*(diag(conj(exp(1j*angle(x)))))';
            L1=L1+trace(inv(Cov_pol)*B{1,mm})*trace(temp1*inv(gamma_MLE));
        end
        pol_detR(ii,jj)=L1+n_pol*log(det(gamma_MLE));
    end
    waitbar(((ii-1)*nc+jj)/(nr*nc),fbar,'in processing¡­');
end

end

