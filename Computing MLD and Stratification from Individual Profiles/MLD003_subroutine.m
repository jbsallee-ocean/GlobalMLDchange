function MLD_subroutine(Source, OceanBasins) 

	disp([Source '---' OceanBasins])
	warning off
	RepData='/net/ether/data/proteo1/jbslod/Taf/BAS/Climatology/Global/Matrix/Update2018/';
	disp('load matrix')

	load([RepData 'Merge_' OceanBasins '_' Source '.mat'])
	disp('end load matrix')
	eval(['Merge=Merge_' OceanBasins '_' Source ';']);
	disp('end load matrix')

	ierror=0;
	disp('create output')
	MLD.length_prof=NaN*ones(size(Merge.lon));
	MLD.pts_above=NaN*ones(size(Merge.lon));
	MLD.pts_below=NaN*ones(size(Merge.lon));
	MLD.gap=NaN*ones(size(Merge.lon));
	MLD.fit=NaN*ones(size(Merge.lon));
	MLD.thrs=NaN*ones(size(Merge.lon));
	MLD.grad=NaN*ones(size(Merge.lon));
	MLD.holte=NaN*ones(size(Merge.lon));
	MLD.perc2a2=NaN*ones(size(Merge.lon));
	MLD.NT15=NaN*ones(size(Merge.lon));
	MLD.NS15=NaN*ones(size(Merge.lon));
	MLD.P=NaN*ones(size(Merge.lon));
	MLD.SA=NaN*ones(size(Merge.lon));
	MLD.CT=NaN*ones(size(Merge.lon));
	MLD.T=NaN*ones(size(Merge.lon));
	MLD.S=NaN*ones(size(Merge.lon));	
	MLD.lon=NaN*ones(size(Merge.lon));
	MLD.lat=NaN*ones(size(Merge.lon));
	MLD.date=NaN*ones(size(Merge.lon));

	sig_grid=20:0.1:30;
	Merge_sig.CT = NaN*ones(length(Merge.lon), length(sig_grid));
	Merge_sig.SA = NaN*ones(length(Merge.lon), length(sig_grid));
	Merge_sig.P = NaN*ones(length(Merge.lon), length(sig_grid));
	Merge_sig.Sigma0=sig_grid;

	disp('start loop')
	for i0=1:length(Merge.T);
		PercentageLoop(i0,length(Merge.T));
		lon_c=Merge.lon(i0); 
		lat_c=Merge.lat(i0); 
		if lon_c>=-180 & lon_c<=360 & lat_c>=-90 & lat_c<=90 
			MLD.lon(i0)=Merge.lon(i0);
			MLD.lat(i0)=Merge.lat(i0);
			date_c=Merge.date(i0); MLD.date(i0)=Merge.date(i0);
			Z_c=Merge.Z;  
			temp=Merge.T(i0,:)'; 
			sal=Merge.S(i0,:)'; 
			inan=isnan(temp+sal);   
			temp(inan)=[];   
			sal(inan)=[];  
			Z_c(inan)=[];    
			pres=gsw_p_from_z(-Z_c,lat_c.*ones(size(Z_c)))';	
			mldindex=1;
			yesplot=0; 
			floatnumber=999;
			MLD.length_prof(i0)=length(temp);	% nbr pts par profils
			[SA, in_ocean]= gsw_SA_from_SP(sal,pres,lon_c,lat_c);
			CT= gsw_CT_from_t(SA,temp,pres);
			density=gsw_sigma0(SA,CT);sigma=density;
	
			test=0
			if test==0
		       	 	ikeep=find(pres<1500);
				if (min(density(ikeep))-max(density(ikeep)))<0.05
					temp=temp(ikeep);
					sal=sal(ikeep);
					pres=pres(ikeep);
					density=density(ikeep);
					sigma=sigma(ikeep);
					CT=CT(ikeep);
					SA=SA(ikeep);
					Z_c=Z_c(ikeep);
				end

	       			ikeep=find(pres<1000);
				if (min(density(ikeep))-max(density(ikeep)))<0.1
					temp=temp(ikeep);
					sal=sal(ikeep);
					pres=pres(ikeep);
					density=density(ikeep);
					sigma=sigma(ikeep);
					CT=CT(ikeep);
					SA=SA(ikeep);
					Z_c=Z_c(ikeep);
				end
			end

			temp_tmp=temp;
			sal_tmp=sal;
			pres_tmp=pres;
			sigma_tmp=sigma;
			if length(temp)>3 & min(pres)<20;
				upperddmax=NaN; 
				mldepthdensp=NaN;
				gdmldp=NaN;
				error_c=NaN;
				try 
					findmld;
					temp=temp_tmp;
					sal=sal_tmp;
					pres=pres_tmp;
					sigma=sigma_tmp;
				catch 
		    			ierror=ierror+1;
		    			error_c=1;
				end

				if isnan(error_c)
					MLD.fit(i0)=upperddmax;
					MLD.thrs(i0)=mldepthdensp;
					MLD.grad(i0)=gdmldp;
					MLD.holte(i0)=mixeddp;
					MLD.T(i0)=mldepthdens_ta;	% temperature moyenne dans la MLD
					MLD.S(i0)=mldepthdens_sa;	% salinité moyenne ds la MLD
					iabove_mld=find(Z_c<mldepthdensp);	% trouve nombre de pts au dessus de MLD de Holte
					ibelow_mld=find(Z_c>=mldepthdensp);	% trouve nombre de pts en dessous de MLD de Holte
					MLD.pts_above(i0)= length(iabove_mld);	% stock nbr de pts au dessus de MLD Holte
					MLD.pts_below(i0)= length(ibelow_mld);	% stock nbr de pts en dessous de MLD Holte
					if ~isempty(ibelow_mld) & ~isempty(iabove_mld);			
						MLD.gap(i0)=abs(min(Z_c(ibelow_mld))-max(Z_c(iabove_mld)));
					else 
						MLD.gap(i0)=NaN;						
					end

					ioksig=find(~isnan(sigma(:)+pres(:)));
					if length(ioksig)>2
						Merge_sig.P(i0,:)=interp1(sigma(ioksig),pres(ioksig),sig_grid);
					end
					ioksig=find(~isnan(sigma(:)+CT(:)));
					if length(ioksig)>2
						Merge_sig.CT(i0,:)=interp1(sigma(ioksig),CT(ioksig),sig_grid);
					end
					ioksig=find(~isnan(sigma(:)+SA(:)));
					if length(ioksig)>2
						Merge_sig.SA(i0,:)=interp1(sigma(ioksig),SA(ioksig),sig_grid);
					end
					CT_20=NaN; CT_200=NaN; CT_base=NaN; CT_baseplus15m=NaN;
					ioksig=find(~isnan(Z_c(:)+CT(:)));
					if length(ioksig)>2
						CT_20=interp1(Z_c(ioksig),CT(ioksig),20);
						CT_200=interp1(Z_c(ioksig),CT(ioksig),200);
						CT_base=interp1(Z_c(ioksig),CT(ioksig),MLD.holte(i0));
						CT_baseplus15m=interp1(Z_c(ioksig),CT(ioksig),MLD.holte(i0)+15);
					end
					SA_20=NaN; SA_200=NaN; SA_base=NaN; SA_baseplus15m=NaN;
					ioksig=find(~isnan(Z_c(:)+SA(:)));
					if length(ioksig)>2				
						SA_20=interp1(Z_c(ioksig),SA(ioksig),20);
						SA_200=interp1(Z_c(ioksig),SA(ioksig),200);
						SA_base=interp1(Z_c(ioksig),SA(ioksig),MLD.thrs(i0));
						SA_baseplus15m=interp1(Z_c(ioksig),SA(ioksig),MLD.thrs(i0)+15);
					end
					MLD.P(i0)=gsw_p_from_z(-abs(MLD.thrs(i0)),lat_c)';
					MLD.SA(i0) = gsw_SA_from_SP(MLD.S(i0),MLD.P(i0),lon_c,lat_c);
					MLD.CT(i0) = gsw_CT_from_t(MLD.SA(i0),MLD.T(i0),MLD.P(i0));

					alpha=gsw_alpha(SA_base,CT_base,MLD.P(i0));
					beta=gsw_beta(SA_base,CT_base,MLD.P(i0));
					MLD.NS15(i0)=beta*((SA_baseplus15m-SA_base)/15);
					MLD.NT15(i0)=alpha*((CT_baseplus15m-CT_base)/15);
					MLD.NS200(i0)=beta*((SA_200-SA_20)/180);
					MLD.NT200(i0)=alpha*((CT_200-CT_20)/180);

					vec=[MLD.fit(i0) MLD.thrs(i0) MLD.holte(i0) ];Perc_2a2(1)=std(vec)./MLD.holte(i0);
					vec=[MLD.grad(i0) MLD.thrs(i0) MLD.holte(i0)];Perc_2a2(2)=std(vec)./MLD.holte(i0);
					vec=[MLD.fit(i0) MLD.grad(i0) MLD.holte(i0) ];Perc_2a2(3)=std(vec)./MLD.holte(i0);
					MLD.perc2a2(i0)=nanmin(Perc_2a2);
				end
			end
		end
	end

	MLD.DateOfCreation=datestr(now);
	eval(['save -v7.3 /net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Database/Update2018/MLD003_' OceanBasins '_' Source '.mat MLD ierror'])
	eval(['save -v7.3 /net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Database/Update2018/MergeSig_' OceanBasins '_' Source '.mat Merge_sig'])

	disp('end this basin, this source')



