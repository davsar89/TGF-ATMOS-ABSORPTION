clear all
close all
clc

set(groot,'defaultfigureposition',[400 250 900 600])
set(0,'DefaultFigureWindowStyle','normal')

MIN_ENERGY = 300; %keV

! rm -rf fused.txt
! cat ./build/output_ascii/* > fused.txt

yy=importdata('fused.txt');

energies = yy(:,6);
source_alts = yy(:,2);
event_nbs = yy(:,4);

ALT_LIST = unique([source_alts; 6; 7; 8]);

% vals = linspace(0,1,64);
% bins = quantile(energies,vals);
bins = logspace(log10(50),log10(40000),32);
bins = sort([bins 511-2 511+2]);

bins=unique(bins);
centers = (bins(1:end-1)+bins(2:end))/2.0;

for ii=1:length(ALT_LIST)
    
    tk = source_alts == ALT_LIST(ii) & energies > MIN_ENERGY;
    
    e_kept = energies(tk);

    if ~isempty( max(event_nbs(tk)))
       nb_sampled(ii) = max(event_nbs(tk));
    else
       nb_sampled(ii) = 0; 
    end
    
    nb_recorded(ii) = sum(tk);

    transmitance(ii) = nb_recorded(ii)/nb_sampled(ii);
    
end

transmitance(isnan(transmitance))=0;

plot(ALT_LIST,transmitance)

set(gca,'yscale','log')

grid on

%title({'TGF photon transmitance (between 0 and 1)', 'as function of photon source altitude.','assuming photon spectrum 1/E*exp(-E/7.3MeV) and Gaussian beaming (sigma 10 deg)','recorded at 500 km altitude, lat = 22deg, long = -77deg','Source photon min energy = 500 keV','Record photon min energy = 300 keV','beaming upwards gaussian sigma 10 deg'})

xlabel('source altitude (km)')
ylabel('transmitance (between 0 and 1)')

saveas(gcf,'result_transmitance_reverse_beam.png')

% text(10,0.015,num2str(transmitance))
%     [N,bins] = histcounts(e_kept,bins);
%     spec{ii} = N./diff(bins);
%     spec{ii} = spec{ii}./trapz(centers,spec{ii});
%     
%     histogram('BinEdges',bins,'BinCounts',spec{ii},'DisplayStyle','stairs','linewidth',2,...
%         'displayname',['E0: ' num2str(ALT_LIST(ii)) ' MeV; ' 'ALPHA: ' num2str(ALPHA_list(ii))])
%     hold on
%     
%     set(gca, 'XScale', 'log')
%     set(gca, 'YScale', 'log')
%     xlabel('energy (keV)')
%     ylabel('n/de spectrum (a.u.)')