% run_megapressproc_GEauto.m
% Jamie Near, McGill University 2017.
%
% USAGE:
% [diffSpec,sumSpec,subSpec1,subSpec2,outw]=run_megapressproc_GEauto(filestring,coilcombos,avgAlignDomain,alignSS);
%
% DESCRIPTION:
% Fully automated processing script for GE MEGA-PRESS MRS data in P-file
% format (GE raw data).  Includes combination of reciever channels,
% removal of bad averages, freqeuncy drift correction, manual alignment of
% edit-on and edit-off spectra, and leftshifting.  This pipeline requires
% no user interaction. This function automatically
% generates an html report to describe the results of each processing step.
%
% INPUTS:
% filestring         = String variable for the name of the P file.
% coilcombos         = (Optional).  A structure obtained by running the
%                      op_getcoilcombos function.  This allows the user to
%                      specify the coil phases and amplitudes as an input,
%                      rather calculating these from the input data by default.
% avgAlignDomain     = (Optional) Perform the spectral registration (drift correction) using
%                     the full spectrum ('t'), or only a limited frequency
%                     range ('f').  Default is 'f'.
% alignSS            = (Optional)
%                     0 - Do not align the edit-on and edit-off subspectra (default).
%                     2 - Perform manual alignment of edit-on and edit-off subspectra.
%
% OUTPUTS:
% diffSpec          = Fully processed difference spectrum.
% sumSpec           = Fully processed sum spectrum.
% subSpec1          = Fully processed MEGA-PRESS subspectrum #1.
% subSpec2          = Fully processed MEGA-PRESS subspectrum #2.
% outw              = Fully processed water unsuppressed spectrum.
%
%function [diffSpec,sumSpec,subSpec1,subSpec2,outw]=run_megapressproc_GEauto(filestring,coilcombos,avgAlignDomain,alignSS);

function preproc_multi_megapress(prefix_path, mat_raww, varargin)


mat_ws = varargin;
processed_specs = cell(1, length(mat_ws));

% process each run up to alignment
for dsi=1:length(mat_ws)
    processed_specs{dsi} = preproc_megapress([prefix_path '_run-' num2str(dsi)], mat_raww, mat_ws{dsi});
end

concat = struct();
for dsi=1:length(mat_ws)
    concat = op_concatAverages(concat, processed_specs{dsi});
end

mkdir(prefix_path);
cd(prefix_path)
%make a new directory for the output report and figures:
mkdir(['./report']);
mkdir(['./report/figs']);

water = 0;
avgAlignDomain = 'f'
out_rm = concat;
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
        iter=0;
        iterin=20;
        p=100;
        fscum=zeros(out_rm.sz(2:end));
        phscum=zeros(out_rm.sz(2:end));
        while (abs(p(1))>0.0003 && iter<iterin)
            iter=iter+1
            close all
            tmax=0.25+0.03*randn(1);
            ppmmin=1.6+0.1*randn(1);
            ppmmaxarray=[3.5+0.1*randn(1,2),4+0.1*randn(1,3),5.5+0.1*randn(1,1)];
            ppmmax=ppmmaxarray(randi(6,1));
            switch avgAlignDomain
                case 't'
                    [out_aa,fs,phs]=op_alignAverages(out_rm2,tmax,'y');
                case 'f'
                    [out_aa,fs,phs]=op_alignAverages_fd(out_rm2,ppmmin,ppmmax,tmax,'y');
                otherwise
                    error('ERROR: avgAlignDomain not recognized!');
            end
            x=repmat([1:size(fs,1)]',1,out_aa.sz(out_aa.dims.subSpecs));
            p=polyfit(x,fs,1)
            fscum=fscum+fs;
            phscum=phscum+phs;
            if driftCorr=='y' || driftCorr=='Y'
                out_rm2=out_aa;
            end
        end
        h=figure('visible','off');
        subplot(2,2,1);
        plot(out_rm.ppm,real(out_rm.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON Before','FontSize',12);
        box off;
        subplot(2,2,2);
        plot(out_aa.ppm,real(out_aa.specs(:,:,1)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-ON After','FontSize',12);
        box off;
        subplot(2,2,3);
        plot(out_rm.ppm,real(out_rm.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF Before','FontSize',12);
        box off;
        subplot(2,2,4);
        plot(out_aa.ppm,real(out_aa.specs(:,:,2)));xlim([1 5]);
        set(gca,'FontSize',8);
        set(gca,'XDir','reverse');
        xlabel('Frequency (ppm)','FontSize',10);
        ylabel('Amplitude(a.u.)','FontSize',10);
        title('Edit-OFF After','FontSize',12);
        box off;
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 20 15]);
        saveas(h,['./report/figs/alignAvgs_prePostFig'],'jpg');
        saveas(h,['./report/figs/alignAvgs_prePostFig'],'fig');
        close(h);
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],fscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Frequency Drift [Hz]','FontSize',10);
        box off;
        legend('Edit-on scans','Edit-off scans','Location','SouthEast');
        legend boxoff;
        title('Estimated Freqeuncy Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,['./report/figs/freqDriftFig'],'jpg');
        saveas(h,['./report/figs/freqDriftFig'],'fig');
        close(h);
        h=figure('visible','off');
        plot([1:out_aa.sz(out_aa.dims.averages)],phscum,'.-','LineWidth',2);
        set(gca,'FontSize',8);
        xlabel('Scan Number','FontSize',10);
        ylabel('Phase Drift [Deg.]','FontSize',10);
        box off;
        legend('Edit-on scans','Edit-off scans');
        legend boxoff;
        title('Estimated Phase Drift','FontSize',12);
        set(h,'PaperUnits','centimeters');
        set(h,'PaperPosition',[0 0 10 10]);
        saveas(h,['./report/figs/phaseDriftFig'],'jpg');
        saveas(h,['./report/figs/phaseDriftFig'],'fig');
        close(h);
        sat='y';
        if sat=='n'
            iter=0;
            p=100;
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
    aligned_spec = out_aa;
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
%now do automatic zero-order phase correction (Use NAA):
out_ls_ss2=op_takesubspec(out_ls,2);
%SpecTool(out_ls_ss2,0.2,1,3.5);
[out_ls_ss2_ph,ph0]=op_autophase(out_ls_ss2,1.9,2.2);
%ph0=input('Input zero-order phase:');
out_ph=op_addphase(out_ls,ph0);
if water
    outw_ph=op_addphase(outw_ls,ph0);
end
%do same phase corection on unprocessed data
out_spectra=out_ph;
on_spectra=op_takesubspec(out_ph,1);
off_spectra=op_takesubspec(out_ph,2);

on_comb=on_spectra;
off_comb=off_spectra;

out_ph = op_concatSubspecs(on_comb, off_comb);
[out_ph,fs,phs]=op_alignAverages(out_ph,.25,'y');
%Now align subspecs if desired:

        out=op_alignMPSubspecs(out_ph);

%aligned_spec = out;
%Make fully processed data;
diffSpec=op_combinesubspecs(out,'diff');
sumSpec=op_combinesubspecs(out,'summ');
subSpec1=op_takesubspec(out,1);
subSpec2=op_takesubspec(out,2);
subSpec2=op_addphase(subSpec2,180);
%Frequency shift all spectra so that Creatine appears at 3.027 ppm:
[subSpec1,frqShift]=op_ppmref(subSpec1,2.9,3.1,3.027);
diffSpec=op_freqshift(diffSpec,frqShift);
sumSpec=op_freqshift(sumSpec,frqShift);
subSpec2=op_freqshift(subSpec2,frqShift);
off_spectra=subSpec1;
on_spectra=subSpec2;
%Make final water unsuppressed data
if water
    if ~isempty(findstr(outw_ph.seq,'edit_529')) || ~isempty(findstr(outw_ph.seq,'jn_svs_special'))
        if outw_ph.dims.subSpecs
            outw=op_combinesubspecs(outw_ph,'diff');
        else
            outw=outw_ph;
        end
    else
        if outw_ph.dims.subSpecs
            outw=op_combinesubspecs(outw_ph,'summ');
        else
            outw=outw_ph;
        end
    end
    outw=op_addphase(outw,-phase(outw.fids(1))*180/pi,0,4.65,1);
else
    outw=0;
end
%Make figure to show the final spectrum:
h=figure('visible','off');
subplot(1,2,1);
plot(subSpec1.ppm,subSpec1.specs,subSpec2.ppm,subSpec2.specs,'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
legend('Edit-ON','Edit-OFF');
legend boxoff;
box off;
title('Result: Subspecs','FontSize',12);
subplot(1,2,2);
plot(diffSpec.ppm,diffSpec.specs,'linewidth',2);xlim([0.2 5.2]);
set(gca,'FontSize',8);
set(gca,'XDir','reverse');
xlabel('Frequency (ppm)','FontSize',10);
ylabel('Amplitude(a.u.)','FontSize',10);
legend('diff');
legend boxoff;
box off;
title('Result: Diff Spectrum','FontSize',12);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition',[0 0 20 10]);
saveas(h,['./report/figs/finalSpecFig'],'jpg');
saveas(h,['./report/figs/finalSpecFig'],'fig');
RF=io_writelcm(diffSpec,['./ave_diff_lcm'],diffSpec.te);


io_writejmrui(diffSpec, sprintf('spectra_jmrui_diff_ave.txt'));
io_writejmrui(subSpec1, sprintf('spectra_jmrui_on_ave.txt'));
io_writejmrui(subSpec2, sprintf('spectra_jmrui_off_ave.txt'));

close all;
