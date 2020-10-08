
for iFit=1:1
%iFit=input('Which Fit?')
iFit
RepFig='/net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/FigureGlobal_revisionFinal2/'
MLDmean                 =0;
MLDtrend                =0;
N2mean                  =0;
N2meanContrib           =0;
N2trend                 =0;
N2trendContrib          =0;
N2_200mean              =0;  
N2_200t                 =0;
Ttrend                  =0;
N2_zonal                =0;
N2_200_zonal            =0;
MLD_zonal               =0;
SSTtrend                =0;
Trend_error             =0;
nomore=1

if nomore==1
load /net/ether/data/proteo1/jbslod/Taf/LOCEAN/MLD/Database/Update2018/MLD003_inventory.mat
NbTotal=squeeze(nansum(squeeze(nansum(Nbdata,1)),1));

%var={'ML_NUM_PROFS',...
    'ML_PRESSURE_WEIGHTED_MEAN','ML_PRESSURE_WEIGHTED_STD','ML_PRESSURE_WEIGHTED_Trend','ML_PRESSURE_WEIGHTED_SE',...
    'ML_N215_WEIGHTED_MEAN','ML_N215_WEIGHTED_STD','ML_N215_WEIGHTED_TREND','ML_N215_WEIGHTED_SE',...
    'ML_NS200_WEIGHTED_MEAN','ML_NS200_WEIGHTED_STD','ML_NS200_WEIGHTED_TREND','ML_NS200_WEIGHTED_SE',...
    'ML_NT200_WEIGHTED_MEAN','ML_NT200_WEIGHTED_STD','ML_NT200_WEIGHTED_TREND','ML_NT200_WEIGHTED_SE',...
    'ML_NS15_WEIGHTED_MEAN','ML_NS15_WEIGHTED_STD','ML_NS15_WEIGHTED_TREND','ML_NS15_WEIGHTED_SE',...
    'ML_NT15_WEIGHTED_MEAN','ML_NT15_WEIGHTED_STD','ML_NT15_WEIGHTED_TREND','ML_NT15_WEIGHTED_SE',...
	'ML_TEMPERATURE_WEIGHTED_MEAN','ML_TEMPERATURE_WEIGHTED_STD','ML_TEMPERATURE_WEIGHTED_Trend','ML_TEMPERATURE_WEIGHTED_SE',...
	'ML_SALINITY_WEIGHTED_MEAN','ML_SALINITY_WEIGHTED_STD','ML_SALINITY_WEIGHTED_Trend','ML_SALINITY_WEIGHTED_SE',...
    'LONGITUDE','LATITUDE','thetadOpt','thetatOpt','thetasOpt','sigmaOpt'};    
%out=Load_File_revision_Final2_newstde(var,iFit)

load MLD_Stratification_1970_2018.mat

NbTotalg=griddata(longrid,latgrid,NbTotal',out.LONGITUDE',out.LATITUDE);

NUMtresh=200;
out.ML_PRESSURE_WEIGHTED_Trend=-out.ML_PRESSURE_WEIGHTED_Trend;
out.ML_PRESSURE_WEIGHTED_MEAN(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_PRESSURE_WEIGHTED_Trend(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_PRESSURE_WEIGHTED_STD(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_PRESSURE_WEIGHTED_SE(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_N215_WEIGHTED_MEAN(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_N215_WEIGHTED_TREND(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_N215_WEIGHTED_STD(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_N215_WEIGHTED_SE(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NS15_WEIGHTED_MEAN(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NS15_WEIGHTED_TREND(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NS15_WEIGHTED_STD(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NS15_WEIGHTED_SE(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NT15_WEIGHTED_MEAN(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NT15_WEIGHTED_TREND(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NT15_WEIGHTED_STD(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NT15_WEIGHTED_SE(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_TEMPERATURE_WEIGHTED_MEAN(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_TEMPERATURE_WEIGHTED_Trend(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_TEMPERATURE_WEIGHTED_STD(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_TEMPERATURE_WEIGHTED_SE(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_SALINITY_WEIGHTED_MEAN(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_SALINITY_WEIGHTED_Trend(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_SALINITY_WEIGHTED_STD(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_SALINITY_WEIGHTED_SE(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NS200_WEIGHTED_MEAN(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NS200_WEIGHTED_TREND(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NS200_WEIGHTED_STD(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NS200_WEIGHTED_SE(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NT200_WEIGHTED_MEAN(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NT200_WEIGHTED_TREND(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NT200_WEIGHTED_STD(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_NT200_WEIGHTED_SE(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_N2200_WEIGHTED_MEAN(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_N2200_WEIGHTED_TREND(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_N2200_WEIGHTED_STD(find(out.ML_NUM_PROFS<NUMtresh))=NaN;
out.ML_N2200_WEIGHTED_SE(find(out.ML_NUM_PROFS<NUMtresh))=NaN;

MLPt=out.ML_PRESSURE_WEIGHTED_Trend;

N2_200=out.ML_NS200_WEIGHTED_MEAN-out.ML_NT200_WEIGHTED_MEAN;
N2_200Error=sqrt(out.ML_NS200_WEIGHTED_STD.^2+out.ML_NT200_WEIGHTED_STD.^2);
N2_200trend=out.ML_NS200_WEIGHTED_TREND-out.ML_NT200_WEIGHTED_TREND;
N2_200SE=sqrt(out.ML_NS200_WEIGHTED_SE.^2+out.ML_NT200_WEIGHTED_SE.^2);
%N2_200trend(find(abs(N2_200trend)<2*N2_200SE))=NaN;

irm=find(abs(out.LATITUDE)<0);


	[NbTotalgok lonn]=change360_2_180(NbTotalg, double(out.LONGITUDE));
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_PRESSURE_WEIGHTED_MEAN(8:10,:,:)))), double(out.LONGITUDE));
	Winter=A;Summer=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_PRESSURE_WEIGHTED_MEAN(1:3,:,:)))), double(out.LONGITUDE));
	Winter(158:end,:)=A(158:end,:);Summer(1:157,:)=A(1:157,:);
	SummerMLD=Summer; WinterMLD=Winter;

    [A lonn]=change360_2_180(double(squeeze((nanmean(out.ML_PRESSURE_WEIGHTED_SE(8:10,:,:)))))*10, double(out.LONGITUDE));
    Winter_SE=A;Summer_SE=A;
    [A lonn]=change360_2_180(double(squeeze((nanmean(out.ML_PRESSURE_WEIGHTED_SE(1:3,:,:)))))*10, double(out.LONGITUDE));
    Winter_SE(158:end,:)=A(158:end,:);Summer_SE(1:157,:)=A(1:157,:);

    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_PRESSURE_WEIGHTED_Trend(8:10,:,:))))*10, double(out.LONGITUDE));
    Winter=A;Summer=A;
    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_PRESSURE_WEIGHTED_Trend(1:3,:,:))))*10, double(out.LONGITUDE));
    Winter(158:end,:)=A(158:end,:);Summer(1:157,:)=A(1:157,:);
	Winter(irm,:)=NaN; Summer(irm,:)=NaN;
    
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_N215_WEIGHTED_MEAN(8:10,:,:)))), double(out.LONGITUDE));
	WinterN2=A;SummerN2=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_N215_WEIGHTED_MEAN(1:3,:,:)))), double(out.LONGITUDE));
	WinterN2(158:end,:)=A(158:end,:);SummerN2(1:157,:)=A(1:157,:);
    
    [A lonn]=change360_2_180(double(squeeze(nanmean(-out.ML_NT15_WEIGHTED_MEAN(8:10,:,:)))), double(out.LONGITUDE));
	WinterNT=A;SummerNT=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(-out.ML_NT15_WEIGHTED_MEAN(1:3,:,:)))), double(out.LONGITUDE));
	WinterNT(158:end,:)=A(158:end,:);SummerNT(1:157,:)=A(1:157,:);

	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS15_WEIGHTED_MEAN(8:10,:,:)))), double(out.LONGITUDE));
	WinterNS=A;SummerNS=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS15_WEIGHTED_MEAN(1:3,:,:)))), double(out.LONGITUDE));
	WinterNS(158:end,:)=A(158:end,:);SummerNS(1:157,:)=A(1:157,:);

	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_N215_WEIGHTED_TREND(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterN2trend=A;SummerN2trend=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_N215_WEIGHTED_TREND(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterN2trend(158:end,:)=A(158:end,:);SummerN2trend(1:157,:)=A(1:157,:);

    [A lonn]=change360_2_180(double(squeeze((nanmean(out.ML_N215_WEIGHTED_SE(8:10,:,:)))))*10, double(out.LONGITUDE));
	WinterN2SE=A;SummerN2SE=A;
	[A lonn]=change360_2_180(double(squeeze((nanmean(out.ML_N215_WEIGHTED_SE(1:3,:,:)))))*10, double(out.LONGITUDE));
	WinterN2SE(158:end,:)=A(158:end,:);SummerN2SE(1:157,:)=A(1:157,:);

    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS15_WEIGHTED_TREND(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterNStrend=A;SummerNStrend=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS15_WEIGHTED_TREND(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterNStrend(158:end,:)=A(158:end,:);SummerNStrend(1:157,:)=A(1:157,:);

	[A lonn]=change360_2_180(double(squeeze(nanmean(-out.ML_NT15_WEIGHTED_TREND(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterNTtrend=A;SummerNTtrend=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(-out.ML_NT15_WEIGHTED_TREND(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterNTtrend(158:end,:)=A(158:end,:);SummerNTtrend(1:157,:)=A(1:157,:);

    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS15_WEIGHTED_SE(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterNSSE=A;SummerNSSE=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS15_WEIGHTED_SE(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterNSSE(158:end,:)=A(158:end,:);SummerNSSE(1:157,:)=A(1:157,:);

	[A lonn]=change360_2_180(double(squeeze(nanmean(-out.ML_NT15_WEIGHTED_SE(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterNTSE=A;SummerNTSE=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(-out.ML_NT15_WEIGHTED_SE(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterNTSE(158:end,:)=A(158:end,:);SummerNTSE(1:157,:)=A(1:157,:);
    
    [A lonn]=change360_2_180(double(squeeze(nanmean(N2_200(8:10,:,:)))), double(out.LONGITUDE));
	WinterN2_200=A;SummerN2_200=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(N2_200(1:3,:,:)))), double(out.LONGITUDE));
	WinterN2_200(158:end,:)=A(158:end,:);SummerN2_200(1:157,:)=A(1:157,:);

	[A lonn]=change360_2_180(double(squeeze(nanmean(N2_200trend(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterN2_200trend=A;SummerN2_200trend=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(N2_200trend(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterN2_200trend(158:end,:)=A(158:end,:);SummerN2_200trend(1:157,:)=A(1:157,:);

    [A lonn]=change360_2_180(double(squeeze(nanmean(N2_200SE(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterN2_200SE=A;SummerN2_200SE=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(N2_200SE(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterN2_200SE(158:end,:)=A(158:end,:);SummerN2_200SE(1:157,:)=A(1:157,:);

    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS200_WEIGHTED_MEAN(8:10,:,:)))), double(out.LONGITUDE));
	WinterNS_200=A;SummerNS_200=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS200_WEIGHTED_MEAN(1:3,:,:)))), double(out.LONGITUDE));
	WinterNS_200(158:end,:)=A(158:end,:);SummerNS_200(1:157,:)=A(1:157,:);

	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS200_WEIGHTED_TREND(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterNS_200trend=A;SummerNS_200trend=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS200_WEIGHTED_TREND(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterNS_200trend(158:end,:)=A(158:end,:);SummerNS_200trend(1:157,:)=A(1:157,:);

    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS200_WEIGHTED_SE(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterNS_200SE=A;SummerNS_200SE=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NS200_WEIGHTED_SE(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterNS_200SE(158:end,:)=A(158:end,:);SummerNS_200SE(1:157,:)=A(1:157,:);
    
    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NT200_WEIGHTED_MEAN(8:10,:,:)))), double(out.LONGITUDE));
	WinterNT_200=A;SummerNT_200=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NT200_WEIGHTED_MEAN(1:3,:,:)))), double(out.LONGITUDE));
	WinterNT_200(158:end,:)=A(158:end,:);SummerNT_200(1:157,:)=A(1:157,:);

	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NT200_WEIGHTED_TREND(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterNT_200trend=A;SummerNT_200trend=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NT200_WEIGHTED_TREND(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterNT_200trend(158:end,:)=A(158:end,:);SummerNT_200trend(1:157,:)=A(1:157,:);

    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NT200_WEIGHTED_SE(8:10,:,:))))*10, double(out.LONGITUDE));
	WinterNT_200SE=A;SummerNT_200SE=A;
	[A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_NT200_WEIGHTED_SE(1:3,:,:))))*10, double(out.LONGITUDE));
	WinterNT_200SE(158:end,:)=A(158:end,:);SummerNT_200SE(1:157,:)=A(1:157,:);
    
    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_TEMPERATURE_WEIGHTED_Trend(8:10,:,:))))*10, double(out.LONGITUDE));
    Winter_Ttrend=A;Summer_Ttrend=A;
    [A lonn]=change360_2_180(double(squeeze(nanmean(out.ML_TEMPERATURE_WEIGHTED_Trend(1:3,:,:))))*10, double(out.LONGITUDE));
    Winter_Ttrend(158:end,:)=A(158:end,:);Summer_Ttrend(1:157,:)=A(1:157,:);

    [A lonn]=change360_2_180(double(squeeze((nanmean(out.ML_TEMPERATURE_WEIGHTED_SE(8:10,:,:)))))*10, double(out.LONGITUDE));
    Winter_TSE=A;Summer_TSE=A;
    [A lonn]=change360_2_180(double(squeeze((nanmean(out.ML_TEMPERATURE_WEIGHTED_SE(1:3,:,:)))))*10, double(out.LONGITUDE));
    Winter_TSE(158:end,:)=A(158:end,:);Summer_TSE(1:157,:)=A(1:157,:);
    
    [Mask lonn]=change360_2_180(double(squeeze(out.ML_NUM_PROFS(2,:,:))), double(out.LONGITUDE));


    Summer(find(NbTotalgok<2 | abs(Summer)<Summer_SE))=NaN;
    SummerN2trend(find(NbTotalgok<2 | abs(SummerN2trend)<SummerN2SE))=NaN;
    Winter(find(NbTotalgok<2 | abs(Winter)<Winter_SE))=NaN;
    WinterN2trend(find(NbTotalgok<2 | abs(WinterN2trend)<WinterN2SE))=NaN;
    SummerNStrend(find(NbTotalgok<2 | abs(SummerNStrend)<SummerNSSE))=NaN;
    SummerNTtrend(find(NbTotalgok<2 | abs(SummerNTtrend)<SummerNSSE))=NaN;

    SummerN2_200trend(find(NbTotalgok<2 | abs(SummerN2_200trend)<SummerN2_200SE))=NaN;
    WinterN2_200trend(find(NbTotalgok<2 | abs(WinterN2_200trend)<SummerN2_200SE))=NaN;


    
    if MLDmean==1
        figure(10); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),WinterMLD)
        caxis([0 400])
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out')
        eval(['print -dpng ' RepFig '/Figure2f_iFit' num2str(iFit) '.png '])
pause(60)

        figure(11); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),SummerMLD)
        caxis([0 100])
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out')
        eval(['print -dpng ' RepFig 'Figure2e_iFit' num2str(iFit) '.png '])
pause(60)
    end
    
    if MLDtrend==1
       figure(12); clf
       Plot_Global_bicolor(lonn,double(out.LATITUDE),Winter)
       caxis([-60 60])
       map=m_colmap('diverging',256);
       h=colorbar('northoutside');
       colormap(map);    
       set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out'); set(h,'xtick',[-60:20:60])
       eval(['print -dpng ' RepFig 'EDFigure6c_iFit' num2str(iFit) '.png '])
pause(60)

        figure(12); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),Summer)
        caxis([-15 15])
        map=m_colmap('diverging',256);
        h=colorbar('northoutside');
        colormap(map)        ;
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out'); set(h,'xtick',[-15:5:15])
        eval(['print -dpng ' RepFig 'Figure3c_iFit' num2str(iFit) '.png '])

	MDT=1;
	if MDT==1
		mdt=Read_Netcdf_JB('/net/ether/data/proteo1/jbslod/Data/Routine-Matlab/MDT_CNES_CLS18.nc',{'longitude','latitude','mean_dynamic_topography'});
		mdt.mean_dynamic_topography(find(mdt.mean_dynamic_topography<-1000))=NaN;
	    	[A lonn]=change360_2_180(mdt.mean_dynamic_topography, mdt.longitude);mdt.mean_dynamic_topography=A*100;mdt.longitude=lonn;
		hold on 
		m_contour(mdt.longitude,mdt.latitude,mdt.mean_dynamic_topography',-100:10:100,'k')
        	eval(['print -dpng ' RepFig 'Figure3c_iFit' num2str(iFit) '_withMDT.png '])
	end


pause(60)      
    end  
    
    if N2mean==1
        figure(14); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),log10(abs(WinterN2)))
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out')
        caxis([-7 -3]);
        set(h,'xtick',[-7:3],'xticklabel',['10^{-7}'; '10^{-6}'; '10^{-5}'; '10^{-4}';  '10^{-3}'])
        eval(['print -dpng ' RepFig 'Figure2d_iFit' num2str(iFit) '.png '])
pause(60)

        figure(15); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),log10(abs(SummerN2)))
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out')
        caxis([-7 -3]);
        set(h,'xtick',[-7:3],'xticklabel',['10^{-7}'; '10^{-6}'; '10^{-5}'; '10^{-4}';  '10^{-3}'])
        eval(['print -dpng ' RepFig 'Figure2c_iFit' num2str(iFit) '.png '])
pause(60)
    end

    if N2meanContrib==1
        figure(16); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),SummerNS./SummerN2*100)
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out');
        caxis([-100 100]);
        map=m_colmap('diverging',256);
        colormap(map)
        eval(['print -dpng ' RepFig '/Figure4b_iFit' num2str(iFit) '.png '])
pause(60)

        figure(17); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),SummerNT./SummerN2*100)
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out');
        caxis([-100 100]);
        map=m_colmap('diverging',256);
        colormap(map)
        eval(['print -dpng ' RepFig '/Figure4a_iFit' num2str(iFit) '.png '])
pause(60)
    end
    

    
    if N2trend==1
        figure(18); clf
        A=SummerN2trend; A(find(abs(A)<10^(-7)))=sign(SummerN2trend(find(abs(A)<10^(-7)))).*10^(-7);
        Plot_Global_bicolor(lonn,double(out.LATITUDE),sign(SummerN2trend).*(7+log10(abs(A))))
        map=m_colmap('diverging',256);
        h=colorbar('northoutside');
        colormap(map)        ;    
        caxis([-3 3])
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out','xtick',-3:3,'xticklabel',['-10^{-4}';'-10^{-5}';'-10^{-6}';'   0    ';' 10^{-6}';' 10^{-5}';' 10^{-4}']);
        eval(['print -dpng ' RepFig '/Figure3b_iFit' num2str(iFit) '.png '])
pause(60)   
        
        figure(19); clf
        A=WinterN2trend; A(find(abs(A)<10^(-7)))=sign(WinterN2trend(find(abs(A)<10^(-7)))).*10^(-7);
        Plot_Global_bicolor(lonn,double(out.LATITUDE),sign(WinterN2trend).*(7+log10(abs(A))))
        map=m_colmap('diverging',256);
        h=colorbar('northoutside');
        colormap(map)        ;    
        caxis([-3 3])
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out','xtick',-3:3,'xticklabel',['-10^{-4}';'-10^{-5}';'-10^{-6}';'   0    ';' 10^{-6}';' 10^{-5}';' 10^{-4}']);
        eval(['print -dpng ' RepFig '/EDFigure6b_iFit' num2str(iFit) '.png '])
pause(60)
    end
 
    if N2trendContrib==1
    	figure(20); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),SummerNStrend./SummerN2trend*100)
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out');
        caxis([-100 100]);
        map=m_colmap('diverging',256);
        colormap(map)
        eval(['print -dpng ' RepFig '/Figure4d_iFit' num2str(iFit) '.png '])
pause(60)

        figure(21); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),SummerNTtrend./SummerN2trend*100)
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out');
        caxis([-100 100]);
        map=m_colmap('diverging',256);
        colormap(map)
        eval(['print -dpng ' RepFig '/Figure4c_iFit' num2str(iFit) '.png '])
pause(60)
    end
    
    if N2_200mean==1
        figure(22); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),log10(abs(WinterN2_200)))
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out')
        caxis([-7 -3]);
        set(h,'xtick',[-7:3],'xticklabel',['10^{-7}'; '10^{-6}'; '10^{-5}'; '10^{-4}';  '10^{-3}'])
        eval(['print -dpng ' RepFig 'Figure2b_iFit' num2str(iFit) '.png '])
pause(60)

        figure(23); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),log10(abs(SummerN2_200)))
        h=colorbar('northoutside');
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out')
        caxis([-7 -3]);
        set(h,'xtick',[-7:3],'xticklabel',['10^{-7}'; '10^{-6}'; '10^{-5}'; '10^{-4}';  '10^{-3}'])
        eval(['print -dpng ' RepFig 'Figure2a_iFit' num2str(iFit) '.png '])
pause(60)
    end
    
    if N2_200t==1
        figure(24); clf
        A=WinterN2_200trend; A(find(abs(A)<10^(-7)))=sign(WinterN2_200trend(find(abs(A)<10^(-7)))).*10^(-7);
        Plot_Global_bicolor(lonn,double(out.LATITUDE),sign(WinterN2_200trend).*(7+log10(abs(A))))
        map=m_colmap('diverging',256);
        h=colorbar('northoutside');
        colormap(map)        ;    
        caxis([-3 3])
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out','xtick',-3:3,'xticklabel',['-10^{-4}';'-10^{-5}';'-10^{-6}';'   0    ';' 10^{-6}';' 10^{-5}';' 10^{-4}']);
        eval(['print -dpng ' RepFig '/EDFigure6a_iFit' num2str(iFit) '.png '])
pause(60)

        figure; clf
        A=SummerN2_200trend; A(find(abs(A)<10^(-7)))=sign(SummerN2_200trend(find(abs(A)<10^(-7)))).*10^(-7);
        Plot_Global_bicolor(lonn,double(out.LATITUDE),sign(SummerN2_200trend).*(7+log10(abs(A))))
        map=m_colmap('diverging',256);
        h=colorbar('northoutside');
        colormap(map)        ;    
        caxis([-3 3])
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out','xtick',-3:3,'xticklabel',['-10^{-4}';'-10^{-5}';'-10^{-6}';'   0    ';' 10^{-6}';' 10^{-5}';' 10^{-4}']);
        eval(['print -dpng ' RepFig '/Figure3a_iFit' num2str(iFit) '.png '])
pause(60)
    end

    if Ttrend==1
    	Summer_Ttrend(find(NbTotalgok<2 ))=NaN;
        figure(29); clf
        Plot_Global_bicolor(lonn,double(out.LATITUDE),Summer_Ttrend)
        caxis([-1 1])
        map=m_colmap('diverging',256);
        h=colorbar('northoutside');
        colormap(map)        ;    
        set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out'); 
        eval(['print -dpng ' RepFig '/EDFigure8a_iFit' num2str(iFit) '.png '])
pause(60)
    end
    
    if MLD_zonal==1
        figure(36); clf; win=10;
        subplot(1,2,2); hold on 
        P=Summer';
        MedianLat=NaN*ones(length(out.LATITUDE),3);
        Nb=P; Nb(find(~isnan(Nb)))=1;Nb=nansum(Nb); P(:,find(Nb<50))=NaN;
        for i=1:length(MedianLat)
            if Nb(i)>=50 & out.LATITUDE(i)<64 & out.LATITUDE(i)>-70
                A=P(:,nanmax(i-win,1):nanmin(i+win,length(out.LATITUDE)));
                MedianLat(i,1)=nanmedian(A(:));
                MedianLat(i,2)=prctile(A(:),66);
                MedianLat(i,3)=prctile(A(:),33);
            end
        end
        XX=[-90 -90 90 90 -90]; YY=[prctile(P(:),33) prctile(P(:),66) prctile(P(:),66) prctile(P(:),33) prctile(P(:),33)];
        plot(MedianLat(:,1),out.LATITUDE,'k','linewidth',2); xlim([-20 20]); ylim([-90 90])
        plot(MedianLat(:,2),out.LATITUDE,'k'); xlim([-20 20]); ylim([-90 90])
        plot(MedianLat(:,3),out.LATITUDE,'k'); xlim([-20 20]); ylim([-90 90])
        plot([0 0],[-90 90],'k'); grid on; grid minor; ylabel('Latitude'); xlabel('m dec^{-1}')
        set(gca,'Fontsize',12,'xtick',[-20:10:20]); %title('Summer MLD')
        eval(['print -dpng ' RepFig '/Figure3c_inset_iFit' num2str(iFit) '.png '])
pause(60)

        figure(37); clf; win=10;
        subplot(1,2,2); hold on 
        P=Winter';
        MedianLat=NaN*ones(length(out.LATITUDE),3);
        Nb=P; Nb(find(~isnan(Nb)))=1;Nb=nansum(Nb); P(:,find(Nb<50))=NaN;
        for i=1:length(MedianLat)
            if Nb(i)>=50& out.LATITUDE(i)<64 & out.LATITUDE(i)>-70
                A=P(:,nanmax(i-win,1):nanmin(i+win,length(out.LATITUDE)));
                MedianLat(i,1)=nanmedian(A(:));
                MedianLat(i,2)=prctile(A(:),66);
                MedianLat(i,3)=prctile(A(:),33);
            end
        end
        disp(['MLD winter:' num2str(nanmedian(P(:)),3) ' m/dec (' num2str(prctile(P(:),33),3) '--' num2str(prctile(P(:),66),3) ')'])
        %-5.65 m/dec (-8.37---3.9)
        XX=[-90 -90 90 90 -90]; YY=[prctile(P(:),33) prctile(P(:),66) prctile(P(:),66) prctile(P(:),33) prctile(P(:),33)];
        %a=fill(YY,XX,'r','edgecolor','none');alpha(a,0.25);
        plot(MedianLat(:,1),out.LATITUDE,'k','linewidth',2); xlim([-20 20]); ylim([-90 90])
        plot(MedianLat(:,2),out.LATITUDE,'k'); xlim([-20 20]); ylim([-90 90])
        plot(MedianLat(:,3),out.LATITUDE,'k'); xlim([-20 20]); ylim([-90 90])
        plot([0 0],[-90 90],'k'); grid on; grid minor; xlabel('m dec^{-1}') %ylabel('Latitude'); 
        set(gca,'Fontsize',12,'xtick',[-20:10:20]);%title('Winter MLD'); ylabel('Latitude'); 
        eval(['print -dpng ' RepFig '/FigureED6c_inset_iFit' num2str(iFit) '.png '])
pause(60)
    end
        
    if N2_zonal==1
        figure(38); clf; win=10;
        subplot(1,2,2); hold on 
        P=SummerN2trend';
        MedianLat=NaN*ones(length(out.LATITUDE),3);
        Nb=P; Nb(find(~isnan(Nb)))=1;Nb=nansum(Nb); P(:,find(Nb<50))=NaN;
        for i=1:length(MedianLat)
            if Nb(i)>=50& out.LATITUDE(i)<64 & out.LATITUDE(i)>-70
                A=P(:,nanmax(i-win,1):nanmin(i+win,length(out.LATITUDE)));
                MedianLat(i,1)=nanmedian(A(:));
                MedianLat(i,2)=prctile(A(:),66);
                MedianLat(i,3)=prctile(A(:),33);
            end
        end
        XX=[-90 -90 90 90 -90]; YY=[prctile(P(:),33) prctile(P(:),66) prctile(P(:),66) prctile(P(:),33) prctile(P(:),33)];
        plot(MedianLat(:,1)*10^5,out.LATITUDE,'k','linewidth',2); ylim([-90 90]); xlim([-2 2]); 
        plot(MedianLat(:,2)*10^5,out.LATITUDE,'k'); ylim([-90 90]); %xlim([-30 30]); 
        plot(MedianLat(:,3)*10^5,out.LATITUDE,'k'); ylim([-90 90]); % xlim([-30 30]); 
        plot([0 0],[-90 90],'k'); grid on; grid minor; ylabel('Latitude'); xlabel('N^{2} trend (10^{-5} s^{-2} dec^{-1})')
        set(gca,'Fontsize',12,'xtick',[-2:1:2]); 
        eval(['print -dpng ' RepFig '/Figure3b_inset_iFit' num2str(iFit) '.png '])
pause(60)

        figure(39); clf; win=10;
        subplot(1,2,2); hold on 
        P=WinterN2trend';
        MedianLat=NaN*ones(length(out.LATITUDE),3);
        Nb=P; Nb(find(~isnan(Nb)))=1;Nb=nansum(Nb); P(:,find(Nb<50))=NaN;
        for i=1:length(MedianLat)
            if Nb(i)>=50& out.LATITUDE(i)<64 & out.LATITUDE(i)>-70
                A=P(:,nanmax(i-win,1):nanmin(i+win,length(out.LATITUDE)));
                MedianLat(i,1)=nanmedian(A(:));
                MedianLat(i,2)=prctile(A(:),66);
                MedianLat(i,3)=prctile(A(:),33);
            end
        end
        XX=[-90 -90 90 90 -90]; YY=[prctile(P(:),33) prctile(P(:),66) prctile(P(:),66) prctile(P(:),33) prctile(P(:),33)];
        plot(MedianLat(:,1)*10^5,out.LATITUDE,'k','linewidth',2); ylim([-90 90]);%xlim([-30 30]); 
        plot(MedianLat(:,2)*10^5,out.LATITUDE,'k'); ylim([-90 90]);
        plot(MedianLat(:,3)*10^5,out.LATITUDE,'k'); ylim([-90 90]);%xlim([-30 30]); 
        plot([0 0],[-90 90],'k'); grid on; grid minor; xlabel('N^{2} trend (10^{-5} s^{-2} dec^{-1})')
        set(gca,'Fontsize',12,'xtick',[-2:1:2]);ylabel('Latitude'); xlim([-2 2]); 
        eval(['print -dpng ' RepFig '/FigureED6b_inset_iFit' num2str(iFit) '.png '])
pause(60)
    end
    

    if N2_200_zonal==1
        figure(38); clf; win=10;
        subplot(1,2,2); hold on 
        P=SummerN2_200trend';
        MedianLat=NaN*ones(length(out.LATITUDE),3);
        Nb=P; Nb(find(~isnan(Nb)))=1;Nb=nansum(Nb); P(:,find(Nb<50))=NaN;
        for i=1:length(MedianLat)
            Nb(i)
            if Nb(i)>=50& out.LATITUDE(i)<64 & out.LATITUDE(i)>-70
                A=P(:,nanmax(i-win,1):nanmin(i+win,length(out.LATITUDE)));
                MedianLat(i,1)=nanmedian(A(:));
                MedianLat(i,2)=prctile(A(:),66);
                MedianLat(i,3)=prctile(A(:),33);
            end
        end
        XX=[-90 -90 90 90 -90]; YY=[prctile(P(:),33) prctile(P(:),66) prctile(P(:),66) prctile(P(:),33) prctile(P(:),33)];
        plot(MedianLat(:,1)*10^7,out.LATITUDE,'k','linewidth',2); ylim([-90 90]); xlim([-30 30]); 
        plot(MedianLat(:,2)*10^7,out.LATITUDE,'k'); ylim([-90 90]); %xlim([-30 30]); 
        plot(MedianLat(:,3)*10^7,out.LATITUDE,'k'); ylim([-90 90]); % xlim([-30 30]); 
        plot([0 0],[-90 90],'k'); grid on; grid minor; ylabel('Latitude'); xlabel('N^{2} trend (10^{-7} s^{-2} dec^{-1})')
        set(gca,'Fontsize',12,'xtick',[-30:10:30]); 
        eval(['print -dpng ' RepFig '/Figure3a_inset_iFit' num2str(iFit) '.png '])
pause(60)

        figure(39); clf; win=10;
        subplot(1,2,2); hold on 
        P=WinterN2_200trend';
        MedianLat=NaN*ones(length(out.LATITUDE),3);
        Nb=P; Nb(find(~isnan(Nb)))=1;Nb=nansum(Nb); P(:,find(Nb<50))=NaN;
        for i=1:length(MedianLat)
            if Nb(i)>=50& out.LATITUDE(i)<64 & out.LATITUDE(i)>-70
                A=P(:,nanmax(i-win,1):nanmin(i+win,length(out.LATITUDE)));
                MedianLat(i,1)=nanmedian(A(:));
                MedianLat(i,2)=prctile(A(:),66);
                MedianLat(i,3)=prctile(A(:),33);
            end
        end
        XX=[-90 -90 90 90 -90]; YY=[prctile(P(:),33) prctile(P(:),66) prctile(P(:),66) prctile(P(:),33) prctile(P(:),33)];
        plot(MedianLat(:,1)*10^7,out.LATITUDE,'k','linewidth',2); ylim([-90 90]);%xlim([-30 30]); 
        plot(MedianLat(:,2)*10^7,out.LATITUDE,'k'); ylim([-90 90]);xlim([-30 30]); 
        plot(MedianLat(:,3)*10^7,out.LATITUDE,'k'); ylim([-90 90]);%xlim([-30 30]); 
        plot([0 0],[-90 90],'k'); grid on; grid minor; xlabel('N^{2} trend (10^{-7} s^{-2} dec^{-1})')
        set(gca,'Fontsize',12,'xtick',[-30:10:30]);ylabel('Latitude'); 
        eval(['print -dpng ' RepFig '/FigureED6a_inset_iFit' num2str(iFit) '.png '])
 pause(60)
   end
    
    
    if SSTtrend==1
        %Goblal Median HadSSTv4 1970-2018: 0.11 deg/dec (0.068--0.153)
        load /net/ether/data/proteo1/jbslod/Data/Wencours/GHRSSTv2/trendSST.mat       
        GHRSSTv2=SummerSST;lon_GHRSSTv2=lon; lat_GHRSSTv2=lat; 
        GHRSSTv2(find(abs(GHRSSTv2)<2*stdSummerSST))=NaN;
        GHRSSTv2(find(GHRSSTv2>-0.005 &GHRSSTv2<0.005))=NaN;
        load /net/ether/data/proteo1/jbslod/Data/Wencours/HadSSTv4/trendSST.mat
        HadSSTv4=SummerSST';lon_hadSST=longitude; lat_hadSST=latitude; 
        SST=Summer_Ttrend;
        SST(find(abs(SST)>1))=NaN;
        
        figure(43)
        g1 = repmat({'SST'},5,1);
        prc25=prctile(SST(:),25)
        prc75=prctile(SST(:),75)
        a=fill([0 4 4 0 0],[prc25 prc25 prc75 prc75 prc25],[0.7 0.7 0.7],'edgecolor','none'); alpha(a,0.4)
        hold on 
        grp=[ones(size(SST(:)')) 2*ones(size(GHRSSTv2(:)')) 3*ones(size(HadSSTv4(:)')) ];
        boxplot([SST(:)' GHRSSTv2(:)' HadSSTv4(:)'  ],grp,'Labels',{'This study (1970-2018)','GHRSSTv2 (1982-2018)','HadSSTv4 (1970-2018)'},'symbol','')
        ylim([-1 1])
        grid on; grid minor;
        ylabel('Temperature trend (^{\circ}C/dec)')
        set(gca,'Fontsize',12)
        eval(['print -dpng ' RepFig '/FigureED8b' num2str(iFit) '.png '])

    end

    if Trend_error==1  
        if iFit<=2
            ERROR1a= Summer_SE;       
            ERROR1b= SummerN2SE;       
            figure(46); clf
            Plot_Global_bicolor(lonn,double(out.LATITUDE),Summer_SE)
            caxis([0 5])
            h=colorbar('northoutside');
            set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out')
            eval(['print -dpng ' RepFig '/MLD003_trendError_summer_iFit' num2str(iFit) '.png '])
            pause(60)
            figure(47); clf
            Plot_Global_bicolor(lonn,double(out.LATITUDE),log10(SummerN2SE))
            caxis([-7 -5])
            h=colorbar('northoutside');
            set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out','xtick',-7:-5,'xticklabel',[' 10^{-7}';' 10^{-6}';' 10^{-5}'])
            eval(['print -dpng ' RepFig 'N2_trendError_summer_iFit' num2str(iFit) '.png '])
            pause(60)
        else
            figure(46); clf
            Plot_Global_bicolor(lonn,double(out.LATITUDE),Summer_SE-ERROR1a)
            caxis([-2 2])
                 map=m_colmap('diverging',256);
                h=colorbar('northoutside');
                colormap(map)        ;    
            set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out')
            eval(['print -dpng ' RepFig '/MLD003_trendError_summer_iFit' num2str(iFit) '.png '])
            pause(60)
            figure(47); clf
            Plot_Global_bicolor(lonn,double(out.LATITUDE),SummerN2SE-ERROR1b)
            caxis([-1 1]*10^(-6))
                map=m_colmap('diverging',256);
                h=colorbar('northoutside');
                colormap(map)        ;    
            set(h,'pos',get(h,'pos')+[.2 .10 -.4 0],'tickdir','out')%,'xtick',-7:-5,'xticklabel',[' 10^{-7}';' 10^{-6}';' 10^{-5}'])
            eval(['print -dpng ' RepFig 'N2_trendError_summer_iFit' num2str(iFit) '.png '])
            pause(60)
        end
    end
    
end
end
