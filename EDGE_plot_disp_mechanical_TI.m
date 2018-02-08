
function EDGE_plot_disp_mechanical_TI(N,par)
    f=par(3);    %free parameters
    A=par(2);
    omega0=par(1);
    load(['EDGE_N' int2str(N) '_f' num2str(f) '_omega0' num2str(omega0) '_A' num2str(A) '.mat'],'k_vec','vals','vecs');
    %load(['EDGE_N_test' int2str(N) '_f' num2str(f) '_omega0' num2str(omega0) '_A' num2str(A) '.mat'],'k_vec','vals','vecs');
    val=(sqrt(-real(vals)));
    NT = 6*N;
    
    figure(1)
    clf
    box on
    set(gcf, 'Units', 'points') 
    xlim([k_vec(1) k_vec(end)])
    xlabel('k','FontSize',16)
    ylabel('\alpha','FontWeight','bold','FontSize',16)
    set(get(gca,'ylabel'),'rotation',0)
    set(get(gca,'ylabel'), 'Position', [-0.35, 16]);
    set(get(gca,'xlabel'), 'Position', [5, 8.7]);
    set(gca,'xTick',0:pi:2*pi)
    set(gca,'xTickLabel',{'0',' ','2\pi'},'FontSize',16)
    
    for jj=1:NT
       
        hold on
        if jj==119 || jj==120 || jj==121 || jj==122 || jj==239 || jj==240 || jj==241 || jj==242
        plot(k_vec,val(:,jj),'-r','LineWidth',3); 
        else
        plot(k_vec,val(:,jj),'-b','LineWidth',0.1); 
    end
    
    ylim([9 19])
    set(gca,'yTick',[9 14 19])
    set(gca,'yTickLabel',{'9','14','19'},'FontSize',16)
end

