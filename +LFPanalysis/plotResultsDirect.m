function plotResultsDirect(Diff,tsec,fvec,MainDir)

    Diffmat = cellfun(@(x) cat(3,x{:}),Diff,'uni',false);
    FIG = figure;
    CS = brewermap(10,'YlOrRd');
    for l = 2:5
        subplot(4,1,l-1);
        hold on;
        line([0 626],[0 0],'linestyle','--','color','k')

        for cnd = 2:2:numel(Diffmat)
            Data = (squeeze(mean(Diffmat{cnd}(:,:,l),1)));
            lg(cnd/2)= plot(tsec,Data,'color',[CS(cnd,:) 1],'linewidth',1.5);

        end
        if l==5
            legend(lg,arrayfun(@(x) ['Cnd-' num2str(x)],2:2:10,'uni',false))
              xlabel('Time(msec)');
            ylabel('Directionality (Upward+)')
        end
        xlim([-100 300])
        ylim([-.1 0.1])
        title(['L' num2str(l)])
        vline(0,'k--');
        set(gca,'fontsize',12)
    end

    set(FIG,'unit','inch','position',[1 1 15 15],'color','w');
    export_fig(FIG,fullfile(MainDir,'Results','Directionality_PDC_allcond'),'-pdf')
    %%
    ind = (tsec>=-100);
    tseccorr = tsec(ind);
    xtickso = find(ismember(tseccorr,-100:100:300));
    ytickso = find(ismember(fvec,20:40:150));


    locs = [1:2:10 2:2:10];
    layers = {2,3,4,5,2:5};
    for l =5:numel(layers)
        FIG2 = figure;
        for cnd = 1:10
            subplot(5,2,locs(cnd));
            imagesc(mean(Diffmat{cnd}(:,ind,layers{l}),3));
            axis xy
            caxis([-.07 .07])
            set(gca,'xtick',xtickso,'xticklabel',tseccorr(xtickso),'ytick',ytickso,'yticklabel',20:40:150);
            vline(find(tseccorr==0),'k--')    
            if cnd==5
               xlabel('time(msec)');
               ylabel('Freq(Hz)')
            end
            title(['Cnd-' num2str(cnd)])
        end
        colormap(jmaColors('coolhotcortex'));
        set(FIG2,'unit','inch','position',[1 1 15 15],'color','w');
    end

    export_fig(FIG2,fullfile(MainDir,'Results','Directionality_PDC_allcond_allfreqs'),'-pdf')
end