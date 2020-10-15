
function outtrend=Sunke_OI_Trend_Final_Sub4_new(yeaTMP , dist_all, B, Omega_tilde,Omega,option);

    %W		= diag(CCC);
    A1 		= [ones(size(yeaTMP))  yeaTMP-2000];
    A2 		= [ones(size(yeaTMP))  yeaTMP-2000 dist_all dist_all.^2];
    
    x=[NaN NaN]; stdx=[NaN NaN];


    if option==1 | option==2
	A=A1;
    elseif option==3 | option==4
	A=A2;
    end


        try 
            %d = eig(Omega_tilde);
            %if all(d > 0
		Omega_tilde_inv = inv(Omega_tilde);
		temp_inv 	= inv(A'*Omega_tilde_inv*A);
		x 		= temp_inv*A'*Omega_tilde_inv*B;
		S 		= temp_inv*A'*Omega_tilde_inv*Omega*Omega_tilde_inv*A*temp_inv;
		stdx 		= sqrt(diag(S));
		x_corr		= sqrt(abs(S(1,2))); %cross correlation of mean (first term) and trend (second term) 
            %else
            %    disp('!!! not positive definite')
            %end
        catch
 	    disp('!!! not positive definite (?)')
            keyboard
        end

    outtrend.mSE=stdx(2); 
    outtrend.mtrend=x(2);
    outtrend.wmeanSE=stdx(1); 
    outtrend.wmean=x(1);
    outtrend.xcorr=x_corr(1);



