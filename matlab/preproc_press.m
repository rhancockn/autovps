% fid_preproc.m for STEAM and PRESS DICOM data
%based on FID-A auto scripts
% Usage:
% fid_preproc(data_ws, data_w, prefix_path)
% raww can be empty
% data_ws
function preproc_press(prefix_path, mat_raww, varargin)

mat_ws = varargin;

iterin=20;
tmaxin=0.2;
aaDomain='f';
water=~isempty(mat_raww);
if ~isempty(mat_raww)
    raww = load(mat_raww);
    raww = raww.svs;
end


spectra=cell(1,length(mat_ws));
for dsi=1:length(mat_ws)

    report_prefix = sprintf('%s/report_%02d', prefix_path,dsi);

    %make a new directory for the output report and figures:
    mkdir(report_prefix );
    mkdir([report_prefix  '/figs']);


    % %read in both datasets:
    raw=load(mat_ws{dsi});
    raw = raw.svs;

    %first step should be to combine coil channels.  To do this find the coil
    %phases from the water unsuppressed data.
    if water
        coilcombos=op_getcoilcombos(raww,1,'h');
        [outw_cc,fidw_pre,specw_pre,phw,sigw]=op_addrcvrs(raww,1,'h',coilcombos);
    else
        coilcombos=op_getcoilcombos(op_averaging(raw),1);
    end
    [out_cc,fid_pre,spec_pre,ph,sig]=op_addrcvrs(raw,1,'h',coilcombos);
    [out_av_cc,fid_av_pre,spec_av_pre]=op_addrcvrs(op_averaging(raw),1,'h',coilcombos);
    raw_av=op_averaging(raw);

    %generate unprocessed spectrum:
    out_noproc=op_averaging(out_cc);
    if water
        outw_noproc=op_averaging(outw_cc);
    end


    %Generate plots showing coil channels before and after phase alignment
    %figure('position',[0 50 560 420]);
    h=figure('visible','off');
    subplot(1,2,1);
    plot(raw_av.ppm,real(raw_av.specs(:,:,1)));xlim([1 5]);
    set(gca,'FontSize',8);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude (a.u.)','FontSize',10);
    title('Before correction','FontSize',12);
    box off;
    subplot(1,2,2);
    plot(raw_av.ppm,real(spec_av_pre(:,:,1)));xlim([1 5]);
    set(gca,'FontSize',12);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    title('After correction','FontSize',12);
    box off;
    set(h,'PaperUnits','centimeters');
    set(h,'PaperPosition',[0 0 20 10]);
    saveas(h,[report_prefix '/figs/coilReconFig'],'jpg');
    saveas(h,[report_prefix '/figs/coilReconFig'],'fig');
    close(h);


    %%%%%%%%OPTIONAL REMOVAL OF BAD AVERAGES FROM FIRST DATASET%%%%%%%%%%%%%%%%%%%%
    close all;
    out_cc2=out_cc;
    nBadAvgTotal=0;
    nbadAverages=1;
    rmbadav='y';
    close all;
    if rmbadav=='n' || rmbadav=='N'
        out_rm=out_cc;
        nsd='N/A';
    else
        sat='n'
        while sat=='n' || sat=='N'
            nsd=4; %Setting the number of standard deviations;
            iter=1;
            nbadAverages=1;
            nBadAvgTotal=0;
            out_cc2=out_cc;
            while nbadAverages>0
                [out_rm,metric{iter},badAverages]=op_rmbadaverages(out_cc2,nsd,'t');
                badAverages;
                nbadAverages=length(badAverages);
                nBadAvgTotal=nBadAvgTotal+nbadAverages;
                out_rm.averages=out_rm.sz(2);
                out_cc2=out_rm;
                out_cc2.averages=out_cc2.sz(2);
                iter=iter+1;
                disp([num2str(nbadAverages) ' bad averages removed on this iteration.']);
                disp([num2str(nBadAvgTotal) ' bad averages removed in total.']);
                close all;
            end
            %figure('position',[0 50 560 420]);
            %Make figure to show pre-post removal of averages
            h=figure('visible','off');
            subplot(1,2,1);
            plot(out_cc.ppm,real(out_cc.specs(:,:)));xlim([1 5]);
            set(gca,'FontSize',8);
            set(gca,'XDir','reverse');
            xlabel('Frequency (ppm)','FontSize',10);
            ylabel('Amplitude(a.u.)','FontSize',10);
            title('Before','FontSize',12);
            box off;
            subplot(1,2,2);
            plot(out_rm.ppm,real(out_rm.specs(:,:)));xlim([1 5]);
            set(gca,'FontSize',8);
            set(gca,'XDir','reverse');
            xlabel('Frequency (ppm)','FontSize',10);
            ylabel('Amplitude(a.u.)','FontSize',10);
            title('After','FontSize',12);
            box off;
            set(h,'PaperUnits','centimeters');
            set(h,'PaperPosition',[0 0 20 15]);
            saveas(h,[report_prefix '/figs/rmBadAvg_prePostFig'],'jpg');
            saveas(h,[report_prefix '/figs/rmBadAvg_prePostFig'],'fig');
            close(h);

            %figure('position',[0 550 560 420]);
            h=figure('visible','off');
            plot([1:length(metric{1})],metric{1},'.r',[1:length(metric{iter-1})],metric{iter-1},'x','MarkerSize',16);
            set(gca,'FontSize',8);
            xlabel('Scan Number','FontSize',10);
            ylabel('Deviation Metric','FontSize',10);
            legend('Before rmBadAv','After rmBadAv');
            legend boxoff;
            title('Deviation Metric','FontSize',12);
            box off;
            set(h,'PaperUnits','centimeters');
            set(h,'PaperPosition',[0 0 20 10]);
            saveas(h,[report_prefix '/figs/rmBadAvg_scatterFig'],'jpg');
            saveas(h,[report_prefix '/figs/rmBadAvg_scatterFig'],'fig');
            close(h);

            %sat1=input('are you satisfied with the removal of bad averages? ','s');
            sat='y';

        end
    end


    %NOW ALIGN AVERAGES:  A.K.A. Frequency Drift Correction.
    driftCorr='y';
    if driftCorr=='n' || driftCorr=='N'
        out_av=op_averaging(out_rm);
        if water
            outw_av=op_averaging(outw_cc);
        end
        fs=0;
        phs=0;
    else
        if water
            outw_aa=op_alignAverages(outw_cc,0.2,'n');
        end
        sat='n';
        out_rm2=out_rm;
        while sat=='n' || sat=='N'
            fsPoly=100;
            phsPoly=1000;
            fscum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
            phscum=zeros(out_rm2.sz(out_rm2.dims.averages),1);
            iter=1;
            while (abs(fsPoly(1))>0.001 || abs(phsPoly(1))>0.01) && iter<iterin
                iter=iter+1
                close all
                tmax=0.25+0.03*randn(1);
                ppmmin=1.6+0.1*randn(1);
                ppmmaxarray=[3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
                ppmmax=ppmmaxarray(randi(6,1));
                switch aaDomain
                    case 't'
                        [out_aa,fs,phs]=op_alignAverages(out_rm2,tmax,'y');
                    case 'f'
                        [out_aa,fs,phs]=op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'y');
                    otherwise
                        error('ERROR: avgAlignDomain not recognized!');
                end

                fsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',fs,1)
                phsPoly=polyfit([1:out_aa.sz(out_aa.dims.averages)]',phs,1)
                iter

                fscum=fscum+fs;
                phscum=phscum+phs;

                if driftCorr=='y' || driftCorr=='Y'
                    out_rm2=out_aa;
                end
            end
            h=figure('visible','off');
            subplot(1,2,1);
            plot(out_rm.ppm,real(out_rm.specs(:,:)));xlim([1 5]);
            set(gca,'FontSize',8);
            set(gca,'XDir','reverse');
            xlabel('Frequency (ppm)','FontSize',10);
            ylabel('Amplitude(a.u.)','FontSize',10);
            title('Before','FontSize',12);
            box off;
            subplot(1,2,2);
            plot(out_aa.ppm,real(out_aa.specs(:,:)));xlim([1 5]);
            set(gca,'FontSize',8);
            set(gca,'XDir','reverse');
            xlabel('Frequency (ppm)','FontSize',10);
            ylabel('Amplitude(a.u.)','FontSize',10);
            title('After','FontSize',12);
            box off;
            set(h,'PaperUnits','centimeters');
            set(h,'PaperPosition',[0 0 20 15]);
            saveas(h,[report_prefix '/figs/alignAvgs_prePostFig'],'jpg');
            saveas(h,[report_prefix '/figs/alignAvgs_prePostFig'],'fig');
            close(h);

            h=figure('visible','off');
            plot([1:out_aa.sz(out_aa.dims.averages)],fscum,'.-','LineWidth',2);
            set(gca,'FontSize',8);
            xlabel('Scan Number','FontSize',10);
            ylabel('Frequency Drift [Hz]','FontSize',10);
            box off;
            legend('Frequency Drift','Location','SouthEast');
            legend boxoff;
            title('Estimated Freqeuncy Drift','FontSize',12);
            set(h,'PaperUnits','centimeters');
            set(h,'PaperPosition',[0 0 10 10]);
            saveas(h,[report_prefix '/figs/freqDriftFig'],'jpg');
            saveas(h,[report_prefix '/figs/freqDriftFig'],'fig');
            close(h);

            h=figure('visible','off');
            plot([1:out_aa.sz(out_aa.dims.averages)],phscum,'.-','LineWidth',2);
            set(gca,'FontSize',8);
            xlabel('Scan Number','FontSize',10);
            ylabel('Phase Drift [Deg.]','FontSize',10);
            box off;
            legend('Phase Drift','Location','SouthEast');
            legend boxoff;
            title('Estimated Phase Drift','FontSize',12);
            set(h,'PaperUnits','centimeters');
            set(h,'PaperPosition',[0 0 10 10]);
            saveas(h,[report_prefix '/figs/phaseDriftFig'],'jpg');
            saveas(h,[report_prefix '/figs/phaseDriftFig'],'fig');
            close(h);

            sat='y';
            if sat=='n'
                iter=0;
                p1=100;
                fscum=zeros(out_rm.sz(2:end));
                phscum=zeros(out_rm.sz(2:end));
                fs2cum=zeros(out_cc.sz(2:end));
                phs2cum=zeros(out_cc.sz(2:end));
                out_rm2=out_rm;
                out_cc2=out_cc;
            end
            totalFreqDrift=mean(max(fscum)-min(fscum));
            totalPhaseDrift=mean(max(phscum)-min(phscum));
            close all
        end
        %now combine the averages averages
        out_av=op_averaging(out_aa);
        if water
            outw_av=op_averaging(outw_aa);
        end
    end

    %now leftshift
    out_ls=op_leftshift(out_av,out_av.pointsToLeftshift);
    if water
        outw_ls=op_leftshift(outw_av,outw_av.pointsToLeftshift);
    end

    %now do automatic zero-order phase correction (Use Creatine Peak):
    [out_ph,ph0]=op_autophase(out_ls,2.9,3.1);
    out_ls_zp=op_zeropad(out_ph,16);

    %And now for water unsuppressed data (use water peak):
    if water
        outw_ls_zp=op_zeropad(outw_ls,16);
        indexw=find(abs(outw_ls_zp.specs)==max(abs(outw_ls_zp.specs(outw_ls_zp.ppm>4 & outw_ls_zp.ppm<5.5))));
        ph0w=-phase(outw_ls_zp.specs(indexw))*180/pi;
        outw_ph=op_addphase(outw_ls,ph0w);
        outw_ls_zp=op_addphase(outw_ls_zp,ph0w);
    end

    %do same phase corection on unprocessed data
    out_noproc=op_addphase(op_leftshift(out_noproc,out_noproc.pointsToLeftshift),ph0);
    if water
        outw_noproc=op_addphase(op_leftshift(outw_noproc,outw_noproc.pointsToLeftshift),ph0w);
    end

    %Frequency shift all spectra so that Creatine appears at 3.027 ppm:
    [~,frqShift]=op_ppmref(out_ls_zp,2.9,3.15,3.027);
    out=op_freqshift(out_ph,frqShift);
    %out=op_freqshift(out_ls_zp,frqShift);
    out_noproc=op_freqshift(out_noproc,frqShift);
    %And now for water unsuppressed data (user water peak and set to 4.65 ppm):
    if water
        [~,frqShiftw]=op_ppmref(outw_ls_zp,4,5.5,4.65);
        outw=op_freqshift(outw_ph,frqShiftw);
        outw_noproc=op_freqshift(outw_noproc,frqShiftw);
    end

    %Make figure to show the final spectrum:
    h=figure('visible','off');
    plot(out.ppm,out.specs,'linewidth',2);xlim([0.2 5.2]);
    set(gca,'FontSize',8);
    set(gca,'XDir','reverse');
    xlabel('Frequency (ppm)','FontSize',10);
    ylabel('Amplitude(a.u.)','FontSize',10);
    legend('summed');
    legend boxoff;
    box off;
    title('Result: Final Spectrum','FontSize',12);
    set(h,'PaperUnits','centimeters');
    set(h,'PaperPosition',[0 0 20 10]);
    saveas(h,[report_prefix '/figs/finalSpecFig'],'jpg');
    saveas(h,[report_prefix '/figs/finalSpecFig'],'fig');

    %write
    io_writejmrui(out, sprintf('%s/spectra_jmrui_%02d.txt', prefix_path,dsi));




    close all;

    %write an html report:
    fid=fopen([report_prefix '/report.html'],'w+');
    fprintf(fid,'<!DOCTYPE html>');
    fprintf(fid,'\n<html>');
    logoPath=which('FID-A_LOGO.jpg');
    fprintf(fid,'\n<img src= " %s " width="120" height="120"></body>',logoPath);
    fprintf(fid,'\n<h1>FID-A Processing Report</h1>');
    fprintf(fid,'\n<h2>Processing pipeline applied to %s data using run_pressproc_auto.m</h2>',raw.seq );
    fprintf(fid,'\n<p>DATE: %s </p>',date);
    fprintf(fid,'\n\n<p> </p>');
    fprintf(fid,'\n\n<h2>Results of multi-coil combination:</h2>');
    fprintf(fid,'\n<img src= "figs/coilReconFig.jpg" width="800" height="400"></body>');
    fprintf(fid,'\n\n<p> </p>');
    fprintf(fid,'\n\n<h2>Results of removal of bad averages:</h2>');
    fprintf(fid,'\n<p>Original number of averages: \t%5.6f </p>',raw.sz(raw.dims.averages));
    fprintf(fid,'\n<p>Number of bad Averages removed:  \t%5.6f </p>',nBadAvgTotal);
    fprintf(fid,'\n<p>Number of remaining averages in processed dataset:  \t%5.6f </p>',out_rm.sz(out_rm.dims.averages));
    fprintf(fid,'\n<p>Bad Averages Removal Threshold was:  \t%2.2f </p>',nsd);
    fprintf(fid,'\n<img src= "figs/rmBadAvg_prePostFig.jpg" width="800" height="600"><img src= "figs/rmBadAvg_scatterFig.jpg " width="800" height="400">');
    fprintf(fid,'\n\n<p> </p>');
    fprintf(fid,'\n\n<p> </p>');
    fprintf(fid,'\n\n<h2>Results of spectral registration:</h2>');
    fprintf(fid,'\n<p>Total frequency drift was: \t%5.6f </p>',max(totalFreqDrift));
    fprintf(fid,'\n<p>Total phase drift was: \t%5.6f </p>',max(totalPhaseDrift));
    fprintf(fid,'\n<img src= "figs/alignAvgs_prePostFig.jpg" width="800" height="600">');
    fprintf(fid,'\n\n<p> </p>');
    fprintf(fid,'\n<img src= "figs/freqDriftFig.jpg" width="400" height="400"><img src="figs/phaseDriftFig.jpg" width="400" height="400">');
    fprintf(fid,'\n\n<p> </p>');
    fprintf(fid,'\n\n<h2>Final Result:</h2>');
    fprintf(fid,'\n<img src= "figs/finalSpecFig.jpg" width="800" height="400">');
    fclose(fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    spectra{dsi}=out;
end

%align and average

if(length(spectra)>1)
    aligned=op_alignAllScans(spectra);
    comb_ave=aligned{1};
    for i=2:length(spectra)
        comb_ave=op_concatAverages(comb_ave, aligned{i});
    end
    comb_ave = op_averaging(comb_ave);
else
    comb_ave=spectra{1};
end

io_writejmrui(comb_ave, sprintf('%s/spectra_jmrui_ave.txt', prefix_path));

fp=fopen([prefix_path '/quality.txt'], 'w+');
fprintf(fp, 'lw\tsnr\n');
fprintf(fp, '%f\t%f\n',op_getLW(comb_ave),op_getSNR(comb_ave));
fclose(fp);
close all;

%Make figure to show the final spectrum:
h=figure('visible','off');
plot(comb_ave.ppm,comb_ave.specs,'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
legend('summed');
legend boxoff;
box off;
title('Result: Combined Spectrum','FontSize',12);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,[prefix_path '/finalSpecFig'],'jpg');
saveas(h,[prefix_path '/finalSpecFig'],'fig');
close all
