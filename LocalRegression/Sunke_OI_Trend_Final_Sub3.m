
function outSub=Sunke_OI_Trend_Final_Sub3(INpara,varINlist,tstep);

	for ii=1:length(varINlist)
		eval([ varINlist{ii} '=INpara.' varINlist{ii}  ';']);
	end


dummy_T = NaN([tim 1],'double');
ML_dens_gr= dummy_T;	
ML_pres_gr= dummy_T;
ML_salt_gr= dummy_T;
ML_temp_gr= dummy_T;
ML_year_gr= dummy_T;
ML_salt_gr_wmean= dummy_T;
ML_temp_gr_wmean= dummy_T;
ML_dens_gr_wmean= dummy_T;
ML_pres_gr_wmean= dummy_T;
ML_salt_gr_wstd= dummy_T;
ML_temp_gr_wstd= dummy_T;
ML_dens_gr_wstd= dummy_T;
ML_pres_gr_wstd= dummy_T;
ML_salt_gr_wSE= dummy_T;
ML_temp_gr_wSE= dummy_T;
ML_dens_gr_wSE= dummy_T;
ML_pres_gr_wSE= dummy_T;
ML_salt_gr_wTrend= dummy_T;
ML_temp_gr_wTrend= dummy_T;
ML_dens_gr_wTrend= dummy_T;
ML_pres_gr_wTrend= dummy_T;
ML_salt_gr_wRMSE= dummy_T;
ML_temp_gr_wRMSE= dummy_T;
ML_dens_gr_wRMSE= dummy_T;
ML_pres_gr_wRMSE= dummy_T;
ML_gr_YearPer5= dummy_T;
ML_gr_YearPer50= dummy_T;
ML_gr_YearPer95= dummy_T;
thetasOpt_wmean= dummy_T;
thetadOpt_wmean= dummy_T;
thetatOpt_wmean= dummy_T;
sigmaOpt_wmean= dummy_T;
for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['ML_' field '_gr_wTrend=dummy_T;']); %=======================MODIF
	eval(['ML_' field '_gr_wRMSE=dummy_T;']); %=======================MODIF
	eval(['ML_' field '_gr_wSE=dummy_T;']); %=======================MODIF
end



                                                
                                                i2_dist = i2_dist_tmp2;
                                                ii_dist = ii_dist_tmp2;

                                                % current month is:
                                                month_step= (tstep/tim)*12; % if tim == 12 then this is equal to month_step = tstep
                                                
                                                % compute fraction of year from month and center all data
                                                % in time around this month set to 0
                                                time_dist=( dyrTMP- (month_step-0.5)/12 ...
                                                    -((dyrTMP-(month_step-0.5)/12) >  0.5)...
                                                    +((dyrTMP-(month_step-0.5)/12) < -0.5) )';
                                                
                                                C_all_T2 = abs(time_dist);
                                                
                                                % EXP 2
                                                C_all_T2 = C_all_T2 .* C_all_T2 ./ (time_scale);
                                                
                                                C_all_T2 = C_all' + C_all_T2;
                                                %ErrorWeight=1./MLpresError(i2_dist)';
                                                %C_all_T = ErrorWeight.*exp(-C_all_T2);
                                                C_all_T = exp(-C_all_T2);
                                                C_all_T = C_all_T .* eps999 + eps100;
						

                                                % only continue if atleast 'min_num_profs' profiles still with us
                                                if numel(ii_dist) > min_num_profs
                                                    
                                                    
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                   % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % COMPUTE MIXED LAYED DEPTH
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    % -------------------------------------------------------------------------
                                                    
                                                    [~, kk ] = sort(C_all_T,'descend');
                                                    
                                                    
                                                    i2_dist = i2_dist(kk);
                                                    %ii_dist = ii_dist(kk);
                                                    
                                                    
                                                    % all ready for computation, do mixed layer and underlying isopycnals:
                                                    sum_str='MLdens(i2_dist) + MLpres(i2_dist) +  MLsalt(i2_dist) + MLtemp(i2_dist)';
                                                    for ifield=1:length(AdditionalFields);
                                                    field=AdditionalFields{ifield};
                                                    eval(['sum_str= [sum_str  ''+ ML' field '(i2_dist)''];']);%=======================	MODIF
                                                    end
                        %%%%

                                                    eval(['index_isfin=find(isfinite(' sum_str ') & dyrMAP(kk)>yeartrendmin);']);                                            
                                                    %index_isfin=find(isfinite(MLdens(i2_dist) + MLatan030(i2_dist) + MLatan0200(i2_dist) + MLpres(i2_dist) +  MLsalt(i2_dist) + MLtemp(i2_dist)));% & thirdofyear');%=============MODIF
                                                    
                                                    ll=numel(index_isfin);
                                                    if ll > min_num_profs
                                                        
                                                        % get number of data points
                                                        % grab closest (weightspace) data profile + random, discard others
                                                        if ll>max_num_profs
                                                            if ll > rnd_within
                                                                ll = rnd_within;
                                                            end
                                                            tmp = seed_rndprofs(seed_rndprofs<=ll);
                                                            index_isfin = index_isfin([1:numsortedprofs tmp(1:num_of_rand_profiles)]);
                                                        end
                                                        
                                                        % compute IQR on MLD density
                                                        
                                                        if doIQR==1
                                                            index_goodIQR  = findgoodIQR(MLdens(i2_dist(index_isfin)),iqr_mult);
                                                            index_goodIQR2 = findgoodIQR(MLpres(i2_dist(index_isfin(index_goodIQR))),iqr_mult);
                                                            index_goodIQR  = index_goodIQR(index_goodIQR2);
                                                        else
                                                            index_goodIQR  = 1:length(index_isfin);
                                                        end

                                                        %  this data is good and can be used
                                                        i2_final = i2_dist(index_isfin(index_goodIQR));
                                                        % ii_final = ii_dist(index_isfin(index_goodIQR));
                                                        i3_final = kk(index_isfin(index_goodIQR));
                                                        % find subset of data for computation:
                                                        
                                                        ll = numel(i2_final);
                                                        
                                                        % set temporary variables
                                                        CaTTMP = C_all_T(i3_final);
                                                        denTMP = MLdens(i2_final);
                                                        preTMP = MLpres(i2_final);
                                                        preErrorTMP = MLpresError(i2_final);
                                                        salTMP = MLsalt(i2_final);
                                                        temTMP = MLtemp(i2_final);
                                                        yeaTMP = dyrMAP(i3_final);
                                                        lonTMP2 = lonTMP(i3_final);
                                                        latTMP2 = latTMP(i3_final);
                                                        for ifield=1:length(AdditionalFields);
                                                            field=AdditionalFields{ifield};
                                                            eval([field 'TMP = ML' field '(i2_final);']); %=======================	MODIF
                                                        end
                                                        mi = (CaTTMP * denTMP) / sum(CaTTMP);
                                                        
                                                        
                                                        % take 100 closest or less
                                                        minleng = min([150; ll]);
                                                        % to compute weighted STD of density

                                                        wfkt = sqrt(var( denTMP(1:minleng),CaTTMP(1:minleng))) .* wfkt_multiplier ;
                                                        ML_Drho(tstep) = wfkt;
                                                        wfkt_sq = wfkt*wfkt;
                                                        
                                                        denTMP2 = abs(denTMP - mi);

                                                              thetadOpt_wmean(tstep)=sqrt(Deccorelation.Space(tstep))*110;
                                                            thetatOpt_wmean(tstep)=sqrt(Deccorelation.Time(tstep))*365;
                                                            thetasOpt_wmean(tstep)=Deccorelation.Phi(tstep);
                                                            sigmaOpt_wmean(tstep)=sqrt(Deccorelation.nugget(tstep));
 %%                                                        
                                                        CCC = 1./CaTTMP.*(Deccorelation.Phi(tstep)+Deccorelation.nugget(tstep)+(preErrorTMP').^2);

                                                        %% Remove outliers
                                                        % too far or too
                                                        % large error
                                                        FieldTrend=[{'den'} {'pre'} {'sal'} {'tem'} {'yea'} AdditionalFields ];
                                                        iok=find(CCC>signal_to_noise);
                                                        %keyboard
                                                        nbobs=length(iok);
                                                        if nbobs>10
                                                            i3_final=i3_final(iok);
                                                            i2_final=i2_final(iok);

                                                            CaTTMP = C_all_T(i3_final);
                                                            denTMP = MLdens(i2_final);
                                                            preTMP = MLpres(i2_final);
                                                            preErrorTMP = MLpresError(i2_final);
                                                            salTMP = MLsalt(i2_final);
                                                            temTMP = MLtemp(i2_final);
                                                            yeaTMP = dyrMAP(i3_final);
                                                            lonTMP2 = lonTMP(i3_final);
                                                            latTMP2 = latTMP(i3_final);
                                                            for ifield=1:length(AdditionalFields);
                                                                field=AdditionalFields{ifield};
                                                                eval([field 'TMP = ML' field '(i2_final);']); %=======================	MODIF
                                                            end
                                                            denTMP2 = abs(denTMP - mi);
                                                            %ErrorWeight=(preTMP'./preErrorTMP').^2;
                                                            %CCC = ErrorWeight.*exp(-(C_all_T2(i3_final) + (((denTMP2.*denTMP2)./wfkt_sq)')));
                                                            CCC = 1./CaTTMP.*(Deccorelation.Phi(tstep)+Deccorelation.nugget(tstep)+(preErrorTMP').^2);

                                                            time_scale_Kuusala=Deccorelation.Time(tstep);
                                                            horiz_scale_Kuusala=Deccorelation.Space(tstep);
                                                            if isnan(time_scale_Kuusala) time_scale_Kuusala=time_scale; end
                                                            if isnan(horiz_scale_Kuusala) horiz_scale_Kuusala=horiz_scale; end

                                                            E_all_T = make_dist_matrix_mex_exp_scal_gauss(...
                                                                    real(dist_vec_complexTMP), ...
                                                                    imag(dist_vec_complexTMP), ...
                                                                    horiz_scale_Kuusala) ...
                                                                    + make_time_matrix_mex_exp_scal_gauss(dyrTMP,time_scale_Kuusala);
                                                             E_Z1 = plus_uptri_and_exp_mex_par2_gauss_new(E_all_T(i3_final,i3_final),denTMP,CaTTMP,CCC); 


                                                            dist_all =abs(dist_vec_complexTMP(i3_final));

                                                            % Generalized Least Square:
                                                            mSEi=NaN; mtrendi=NaN; mRMSEi=NaN;mSEs=NaN;mtrends=NaN; mRMSEs=NaN; 
                                                            mSEy=NaN; mtrendy=NaN; mRMSEy=NaN;mSEt=NaN; mtrendt=NaN; mRMSEt=NaN;

                                                            for ifield=1:length(FieldTrend);
                                                                field=FieldTrend{ifield};
                                                                eval(['INtrend{ifield} = ' field 'TMP;']) 
                                                            end
                                                            %tic 
                                                            %keyboard
                                                            optionFit=optionFit;
                                                            %parfor ifield=1:length(FieldTrend);
                                                            for ifield=1:length(FieldTrend);
                                                                    outtrend{ifield}=Sunke_OI_Trend_Final_Sub4(yeaTMP , dist_all, INtrend{ifield}, E_Z1,CCC,optionFit);
                                                            end
                                                            %toc
                                                            clear INtrend							
                                                            for ifield=1:length(FieldTrend);
                                                                if ifield==1	mSEi=outtrend{ifield}.mSE; mtrendi=outtrend{ifield}.mtrend; end	
                                                                if ifield==2	mSEp=outtrend{ifield}.mSE; mtrendp=outtrend{ifield}.mtrend; end	
                                                                if ifield==3	mSEs=outtrend{ifield}.mSE; mtrends=outtrend{ifield}.mtrend; end	
                                                                if ifield==4	mSEt=outtrend{ifield}.mSE; mtrendt=outtrend{ifield}.mtrend; end	
                                                                if ifield==5	mSEy=outtrend{ifield}.mSE; mtrendy=outtrend{ifield}.mtrend; end	

                                                                if ifield==1	mi=outtrend{ifield}.wmean; mstdi=outtrend{ifield}.wmeanSE; end	
                                                                if ifield==2	mp=outtrend{ifield}.wmean; mstdp=outtrend{ifield}.wmeanSE; end	
                                                                if ifield==3	ms=outtrend{ifield}.wmean; mstds=outtrend{ifield}.wmeanSE; end	
                                                                if ifield==4	mt=outtrend{ifield}.wmean; mstdt=outtrend{ifield}.wmeanSE; end	
                                                                if ifield==5	my=outtrend{ifield}.wmean; mstdy=outtrend{ifield}.wmeanSE; end	

                                                                if ifield>5	
                                                                    field=FieldTrend{ifield};
                                                                    eval(['mSEi' field '=outtrend{ifield}.mSE; mtrendi' field '=outtrend{ifield}.mtrend;']);
                                                                    eval(['mi' field '=outtrend{ifield}.wmean; mstdi' field '=outtrend{ifield}.wmeanSE;']);
                                                                end
                                                            end
                                                        else 
                                                            if isempty(find(~isnan(CCC)))
                                                                disp('no data: NaN')
                                                            else
                                                                disp('not enough data')
                                                            end
                                                            
                                                             for ifield=1:length(FieldTrend);
                                                                if ifield==1	mSEi=NaN; mtrendi=NaN; end	
                                                                if ifield==2	mSEp=NaN; mtrendp=NaN; end	
                                                                if ifield==3	mSEs=NaN; mtrends=NaN; end	
                                                                if ifield==4	mSEt=NaN; mtrendt=NaN; end	
                                                                if ifield==5	mSEy=NaN; mtrendy=NaN; end	

                                                                if ifield==1	mi=NaN; mstdi=NaN; end	
                                                                if ifield==2	mp=NaN; mstdp=NaN; end	
                                                                if ifield==3	ms=NaN; mstds=NaN; end	
                                                                if ifield==4	mt=NaN; mstdt=NaN; end	
                                                                if ifield==5	my=NaN; mstdy=NaN; end	

                                                                if ifield>5	
                                                                    field=FieldTrend{ifield};
                                                                    eval(['mSEi' field '=NaN; mtrendi' field '=NaN;']);
                                                                    eval(['mi' field '=NaN; mstdi' field '=NaN;']);
                                                                end
                                                             end
                                                        end


                                                        % remove  mean
                                                        denTMP = (denTMP - mi);
                                                        preTMP = (preTMP - mp);
                                                        salTMP = (salTMP - ms);
                                                        temTMP = (temTMP - mt);
                                                        
                                                        for ifield=1:length(AdditionalFields);
                                                            field=AdditionalFields{ifield};
                                                            eval([field 'TMP = (' field 'TMP - mi' field ');']); %=======================	MODIF
                                                        end

                                                        % remove  mean for year: 
                                                        yeaTMP = (yeaTMP - my);

                                                        % write weighted mean of dens and pres to file
                                                        ML_dens_gr_wmean(tstep) = mi;
                                                        ML_pres_gr_wmean(tstep) = mp;
                                                        ML_salt_gr_wmean(tstep) = ms;
                                                        ML_temp_gr_wmean(tstep) = mt;
                                                        for ifield=1:length(AdditionalFields);
                                                            field=AdditionalFields{ifield};
                                                            eval(['ML_' field '_gr_wmean(tstep) = mi' field ';']); %=======================	MODIF
                                                        end
                                                        
                                                        ML_dens_gr_wstd(tstep) = mstdi;
                                                        ML_pres_gr_wstd(tstep) = mstdp;
                                                        ML_salt_gr_wstd(tstep) = mstds;
                                                        ML_temp_gr_wstd(tstep) = mstdt;
                                                        for ifield=1:length(AdditionalFields);
                                                            field=AdditionalFields{ifield};
                                                            eval(['ML_' field '_gr_wstd(tstep) = mstdi' field ';']); %=======================	MODIF
                                                        end
                                                        
                                                        ML_dens_gr_wSE(tstep) = mSEi;
                                                        ML_pres_gr_wSE(tstep) = mSEp;
                                                        ML_salt_gr_wSE(tstep) = mSEs;
                                                        ML_temp_gr_wSE(tstep) = mSEt;
                                                        for ifield=1:length(AdditionalFields);
                                                            field=AdditionalFields{ifield};
                                                            eval(['ML_' field '_gr_wSE(tstep) = mSEi' field ';']); %=======================	MODIF
                                                        end
                                                        
                                                        ML_dens_gr_wTrend(tstep) = mtrendi;
                                                        ML_pres_gr_wTrend(tstep) = mtrendp;
                                                        ML_salt_gr_wTrend(tstep) = mtrends;
                                                        ML_temp_gr_wTrend(tstep) = mtrendt;
                                                        for ifield=1:length(AdditionalFields);
                                                            field=AdditionalFields{ifield};
                                                            eval(['ML_' field '_gr_wTrend(tstep) = mtrendi' field ';']); %=======================	MODIF
                                                        end
                                                        
                                                        ML_dens_gr_wRMSE(tstep) = NaN;
                                                        ML_pres_gr_wRMSE(tstep) = NaN;
                                                        ML_salt_gr_wRMSE(tstep) = NaN;
                                                        ML_temp_gr_wRMSE(tstep) = NaN;
                                                        for ifield=1:length(AdditionalFields);
                                                            field=AdditionalFields{ifield};
                                                            eval(['ML_' field '_gr_wRMSE(tstep) = NaN;']); %=======================	MODIF
                                                        end
							
                                                        ML_gr_YearPer5(tstep) = prctile(dyrMAP,5);
                                                        ML_gr_YearPer50(tstep) = prctile(dyrMAP,50);
                                                        ML_gr_YearPer95(tstep) = prctile(dyrMAP,95);

                                                        CCC = 1./CaTTMP.*(Deccorelation.Phi(tstep)+Deccorelation.nugget(tstep)+(preErrorTMP').^2);                                                        
                                                        %CCC = CCC .* eps999 + eps100;
                                                        tmpsw = sum(CCC);
                                                        
                                                        % save summed weight to file
                                                        ML_weight(tstep) = tmpsw;
        
                                                        ML_salt_gr(tstep) = NaN;
                                                        ML_temp_gr(tstep) = NaN;
                                                        ML_dens_gr(tstep) = NaN;
                                                        ML_pres_gr(tstep) = NaN;
                                                        ML_year_gr(tstep) = NaN;
                                                        for ifield=1:length(AdditionalFields);
                                                            field=AdditionalFields{ifield};
                                                            eval(['ML_' field '_gr(tstep) = NaN;']); %=======================	MODIF
                                                        end                                                        
                                                        %save all errors into grid point TIME matrix
                                                        ML_area(tstep) = NaN;
                                                        ML_error(tstep) = NaN;

                                                        ML_r_dist(tstep) = (CaTTMP * dist_vec_complexTMP(i3_final)) *110;
                                                        old_depth = ML_pres_gr(tstep);
                                                        
                                                        
                                                        
                                                    else
                                                        old_depth = 0;
                                                    end
                                                    
                                                    
                                                    if isfinite(ML_dens_gr(tstep) + ML_pres_gr(tstep)) && ML_dens_gr(tstep) > 2
                                                        sg_ind_ML = find(sg0_pleth < ML_dens_gr(tstep),1,'last'); % start at this isopycnal IN ML
                                                        sg_ind = sg_ind_ML;
                                                        sgindmax = 1;
                                                    else
                                                        sg_ind_ML = 1;
                                                        sg_ind=1;
                                                        sgindmax = 1;
                                                    end
                                                    

                                                end

varlist={'ML_dens_gr';	'ML_pres_gr';'ML_salt_gr';'ML_temp_gr';'ML_year_gr';'ML_salt_gr_wmean';'ML_temp_gr_wmean';'ML_dens_gr_wmean';'ML_pres_gr_wmean';'ML_salt_gr_wstd';'ML_temp_gr_wstd';'ML_dens_gr_wstd';'ML_pres_gr_wstd';'ML_salt_gr_wSE';'ML_temp_gr_wSE';'ML_dens_gr_wSE';
'ML_pres_gr_wSE';'ML_salt_gr_wTrend';'ML_temp_gr_wTrend';'ML_dens_gr_wTrend';'ML_pres_gr_wTrend';'ML_salt_gr_wRMSE';'ML_temp_gr_wRMSE';'ML_dens_gr_wRMSE';'ML_pres_gr_wRMSE';'ML_gr_YearPer5';'ML_gr_YearPer50';'ML_gr_YearPer95';'thetasOpt_wmean';
'thetadOpt_wmean';'thetatOpt_wmean';'sigmaOpt_wmean';'ML_num_used'};

for ifield=1:length(AdditionalFields);
	field=AdditionalFields{ifield};
	eval(['varlist= [varlist ; {''ML_' field '_gr_wstd''}];']); %=======================MODIF
	eval(['varlist= [varlist ; {''ML_' field '_gr_wmean''}];']); %=======================MODIF
	eval(['varlist= [varlist ; {''ML_' field '_gr_wTrend''}];']); %=======================MODIF
	eval(['varlist= [varlist ; {''ML_' field '_gr_wRMSE''}];']); %=======================MODIF
	eval(['varlist= [varlist ; {''ML_' field '_gr_wSE''}];']); %=======================MODIF
end


for ivar=1:length(varlist)
	eval(['outSub.' varlist{ivar} '=' varlist{ivar} ';']);
end




