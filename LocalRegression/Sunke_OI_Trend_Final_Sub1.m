
function out=Sunke_OI_Trend_Final_Sub1(INpara,varINlist,istep);



for ii=1:length(varINlist)
	eval([ varINlist{ii} '=INpara.' varINlist{ii}  ';']);
end

jstep=1;
        disp([' work on latitude: ', ...
            num2str(y(istep,jstep)), '?N  ...'])
	eval(['!echo work on latitude: ' num2str(y(istep,jstep)) '?N  ... >> ' logfile])           


    tic
    if y(istep,1) >= continue_at && y(istep,1) <= continue_upto
        time_now = rem(now,1);


	INpara=INpara;
	jstep=jstep;
	%parfor jstep=1:n %longitude loop 360 times / 720 times if 0.5? ...	
	for jstep=1:n %longitude loop 360 times / 720 times if 0.5? ...	
		out{jstep}=Sunke_OI_Trend_Final_Sub2(INpara,varINlist,jstep,istep);
	end

	jstep=1;	
        % print out latitude processed and time used since start ...
        fprintf(' \n ')
        disp([' calculated latitude: ', ...
            num2str(y(istep,jstep)), '?N  - running since ',...
            datestr(now-tstart,'HH:MM:SS'), ' HH:MM:SS.'])

	eval(['!echo calculated latitude: ' num2str(y(istep,jstep)) '?N  - running since '  datestr(now-tstart,'HH:MM:SS') ' HH:MM:SS. >> ' logfile])           
    else

	for jstep=1:n %longitude loop 360 times / 720 times if 0.5? ...	
		out{jstep}.istep=NaN;
		out{jstep}.jstep=NaN;
	end

    end
    toc

