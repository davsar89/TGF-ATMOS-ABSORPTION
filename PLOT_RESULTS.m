clear all
close all
clc

E0_list = [7.3 7.6 7.9];
ALPHA_list = [1 1.2 1.2];

! cat /Home/siv29/dsa030/Desktop/GEANT4/CHECK_TGF_SPECTRUM_AFTER_PROPAGATION_ANDERS_PAPER/build/output_ascii/* > fused.txt

yy=importdata('fused.txt');

energies = yy(:,6);
E0 = yy(:,13);
ALPHA = yy(:,14);

bins=logspace(log10(50),log10(40000),128);

for ii=1:length(E0_list)
    
    tk = E0 == E0_list(ii) & ALPHA == ALPHA_list(ii);
    
    e_kept = energies(tk);
    
    [N,bins] = histcounts(e_kept,bins);
    spec{ii} = N./diff(bins);
    spec{ii} = spec{ii}./max(spec{ii});
    
    histogram('BinEdges',bins,'BinCounts',spec{ii},'DisplayStyle','stairs',...
        'displayname',['E0: ' num2str(E0_list(ii)) ' MeV; ' 'ALPHA: ' num2str(ALPHA_list(ii))])
    hold on
    
    set(gca, 'XScale', 'log')
    set(gca, 'YScale', 'log')
    xlabel('energy (keV)')
    ylabel('n/de spectrum (a.u.)')
    
end

legend('show')

grid on

set(gca,'xlim',[45 41000])

title('TGF at 15 km altitude, record at 400 km altitude, angular dist. Gaussian with sigma = 20')

saveas(gcf,'result_compa_usual_dwyer_geant4.png')