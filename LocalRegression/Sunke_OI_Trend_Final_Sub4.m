
function outtrend=Sunke_OI_Trend_Final_Sub4(yeaTMP , dist_all, YY, V1, CCC, option);

    Am1=[ones(size(yeaTMP))];
    Am2=[ones(size(yeaTMP)) dist_all dist_all.^2];
    A1=[ones(size(yeaTMP))  yeaTMP-2000];
    A2=[ones(size(yeaTMP))  yeaTMP-2000 dist_all dist_all.^2];
    x=[NaN NaN]; stdx=[NaN NaN];
    if option==1
        try 
            d = eig(V1);
            if all(d > 0)
                [x,stdx] = lscov(A1,YY,V1); %disp('positive definite')
            elseif all(d >= 0)
                x=10^(10)*[1 1];  disp('semi positive definite')          
            elseif all(d <= 0)
                x=-10^(10)*[1 1]; disp('negative definite')
            else
                disp('indefinite')
            end
        catch
            keyboard
        end
    elseif option==2
        V2 =diag(CCC);    
        [x,stdx] = lscov(A1,YY,V2);
    elseif option==3
        [R,p] = chol(V1);
        if p==0
            [x,stdx] = lscov(A2,YY,V1);
        end
    elseif option==4
        V2 =diag(CCC);    
        [x,stdx] = lscov(A2,YY,V2);
    elseif option==5
        linearModel = fitlm(yeaTMP , YY, 'linear','Weights',CCC);
        x(2)=linearModel.Coefficients{2,1};
        tmpsw = sum(CCC); x(1) = (CCC * YY) /tmpsw;
        stdx(2)=linearModel.Coefficients{2,2}; 
        stdx(1) = sqrt((CCC * (YY - x(1)).^2) /tmpsw);
    elseif option==6
        [x,stdx] = lscov(Am1,YY,V1); %disp('positive definite')
        xtmp=x(1); stdxtmp=stdx(1);
        [x,stdx] = lscov(Am2,YY,V1); %disp('positive definite')
        xtmp2=x(1); stdxtmp2=stdx(1);
        x=[xtmp xtmp2]; stdx=[stdxtmp stdxtmp2];
    end
    outtrend.mSE=stdx(2); 
	outtrend.mtrend=x(2);
	outtrend.wmeanSE=stdx(1); 
	outtrend.wmean=x(1);

%     A=[ones(size(yeaTMP))  (yeaTMP-2000) dist_all dist_all.^2];
%     [x,stdx] = lscov(A,YY,V);
% 	outtrend.mSE=stdx(2); 
% 	outtrend.mtrend=x(2);
% 	outtrend.wmeanSE=stdx(1); 
% 	outtrend.wmean=x(1);
% 
% 
% 
%  	tic; linearModel = fitlm(yeaTMP , YY, 'linear','Weights',CCC);toc
%  	outtrend.mSE=linearModel.Coefficients{2,2}; 
%  	outtrend.mtrend=linearModel.Coefficients{2,1}; 
%  	outtrend.mRMSE=linearModel.RMSE;	
%     toc


