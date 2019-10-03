function plot_PDC_Directionality(Direction_Stats,tsec,fvec,tWin,labels,ROIs,opt,SaveFigName,StatType)

% this function plots directionality results:
% Works with bootstrapping results!

%%
FS  =   12;
load LayerColors.mat 
lnum     =   numel(labels); %%%%%
for roi     =   1:numel(ROIs) % for each hemi
   

   %-------------------------Averaged Over Frequencies-----------------
    Fig2 = figure;
    clear SP;
    line([tWin(1) tWin(2)],[0 0],'linestyle','--','color','k','linewidth',1.3);
    hold on;
    for i = 2:lnum-1
        LIdx = find(strcmp({Direction_Stats.LName},[labels{i} '_' ROIs{roi}]));
        SP(i-1) = plot(tsec,sum(Direction_Stats(LIdx).Median,1),'color',Colors(i,:),'linewidth',2);
    end
    xlim(tWin)
    xlabel('Time (ms)')
    vline([0],'k--')
    legend(SP,labels(2:lnum-1));
    set(gca,'fontsize',16);
    set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
    title ('Upward flow of layers')
    export_fig(Fig2,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AveragedFreqs_' ROIs{roi}]),'-pdf');

    %--------------------------All Frequencies-------------------------
    ind         =       (tsec>=tWin(1)) & (tsec<=tWin(2)); % time window for plotting
    tseccorr    =       tsec(ind);
    Fig22       =       figure;
    tsec2       =       (tsec>=0);
    for i = 2:lnum-1

        SP(i-1)     =       subplot(4,1,i-1);
        LIdx        =       find(strcmp({Direction_Stats.LName},[labels{i} '_' ROIs{roi}]));
        Dataplot    =       Direction_Stats(LIdx).Median(:,ind);
        TP = zeros(size(Direction_Stats(LIdx).Median));
        if strcmp(StatType,'bootstrap')
            Transplot   =       Direction_Stats(LIdx).PostStim.h;          
            TP(:,tsec2) = Transplot;
        elseif strcmp(StatType,'permutation')
            SigCs = [Direction_Stats(LIdx).Clusters.p{:}]<.05; % Find significant clusters
            Nodes   =       Direction_Stats(LIdx).Clusters.nodes(SigCs); 
            Nodes = cat(1,Nodes{:});
            TP(Nodes) = 1;
        end
        TP = TP(:,ind);
        imagesc(Dataplot,'alphadata',TP*.7+.3);
        axis xy;
        interv  = round(round((tWin(2)-tWin(1))/5)/10)*10;
        a = .2;
        par1 = 0:-interv:tWin(1);
        xtickso = find(ismember(round(tseccorr/a)*a,[par1(end:-1:1) 0:interv:tWin(2)]));
        xtickso(diff(xtickso)==1)=[];
        set(gca,'xtick',xtickso,'xticklabel',round(tseccorr(xtickso)),'ytick',1:40:numel(fvec),'yticklabel',fvec(1:40:end));
        M = max(abs(Dataplot(:)));
        caxis([-M M]);
        %caxis([-.2 .2]);
        vline(xtickso(tseccorr(xtickso)==0),'k--');
        colorbar;
        if i==5
            xlabel('Time(msec)');
            ylabel('Frequency(Hz)');
        end
        title(labels{i})
        set(gca,'fontsize',FS)
    end
    colormap(jmaColors('coolhotcortex'))
    set(Fig22,'unit','inch','position',[1 1 10 15],'color','w')

    export_fig(Fig22,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AllFreqs_' ROIs{roi}]),'-pdf');

    %--------------------average over all Layers-----------------------

    LIdx        =       find(strcmp({Direction_Stats.LName},['All_' ROIs{roi}]));
    Dataplot    =       Direction_Stats(LIdx).Median(:,ind);
    TP = zeros(size(Direction_Stats(LIdx).Median));
    
    if strcmp(StatType,'bootstrap')
        Transplot   =       Direction_Stats(LIdx).PostStim.h;          
        TP(:,tsec2) = Transplot;
    elseif strcmp(StatType,'permutation')
        SigCs = [Direction_Stats(LIdx).Clusters.p{:}]<.05; % Find significant clusters
        Nodes   =       Direction_Stats(LIdx).Clusters.nodes(SigCs); 
        Nodes = cat(1,Nodes{:});
        TP(Nodes) = 1;
    end

    TP = TP(:,ind);

    Fig3 = figure;
    imagesc(Dataplot,'alphadata',TP*.7+.3);
    axis xy;
    interv  = round(round((tWin(2)-tWin(1))/5)/10)*10;
    a = .2;
    par1 = 0:-interv:tWin(1);
    xtickso = find(ismember(round(tseccorr/a)*a,[par1(end:-1:1) 0:interv:tWin(2)]));
    xtickso(diff(xtickso)==1)=[];
    set(gca,'xtick',xtickso,'xticklabel',tseccorr(xtickso),'ytick',1:40:numel(fvec),'yticklabel',fvec(1:40:end));
    M = max(abs(Dataplot(:)));
    caxis([-M M]);
    vline(xtickso(tseccorr(xtickso)==0),'k--');
    colorbar;
    xlabel('Time(msec)');
    ylabel('Frequency(Hz)');
    set(gca,'fontsize',FS);
    colormap(jmaColors('coolhotcortex'))
    set(Fig3,'unit','inch','position',[1 1 10 4],'color','w');
    export_fig(Fig3,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AveragedLayers_' ROIs{roi}]),'-pdf');

    close all

end
end