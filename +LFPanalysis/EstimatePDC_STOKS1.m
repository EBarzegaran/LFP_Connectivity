function [PDC, fvec] = EstimatePDC_STOKS1(Projfolder,varargin)

%% set default values
    opt = ParseArgs(varargin, ...
        'plotfig'       ,true, ...
        'ModOrds'       ,15,...
        'Cnd'           ,6, ...
        'recalc'        ,false, ...
        'Normalize'     ,'All',...
        'Freqband'      ,[20 150],...
        'SStat'         ,'Paired',...
        'doStats'       ,false,...
        'NormPrestim'   ,false,...
        'doPermute'     ,false,...
        'typePermute'   ,'layers',...
        'NPerm'         ,100,...
        'TimeWin'       ,[-100 100],...
        'figpath'       ,fullfile(Projfolder,'Results') ...
        );

    animals = dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals = {animals.name};
    animID = cellfun(@(x) x(1:4),animals,'uni',false); % animal IDs
    
    opt.Freqband = round(opt.Freqband); % Round the frequency band for now
    SaveFileName = ['PDCSTOK_MORD' num2str(opt.ModOrds)];
%% Estimate MVAR parameters using STOK algorithm
    xticklabels = -100:100:100;
    xticks = linspace (0,1250,numel(xticklabels));
    
    if ~exist(fullfile(Projfolder,[SaveFileName '.mat']),'file') || opt.recalc
        for subj = 1:numel(animals)
            
            %---------------------- Load the Data --------------------
            disp(animID{subj});
            LFP = load(fullfile(Projfolder,animals{subj}),'lfpRat');
            load(fullfile(Projfolder,animals{subj}),'tsec');
           
            %--------------- Estimate MVAR params and PDC -----------------
            % (1) Prepare the data
            epochs = permute(LFP.lfpRat(:,:,:), [3,2,1]); %trials, nodes, time
            % (2) prepare the parameters
            ff=.99;
            keepdiag = 1; % 
            measure='PDCnn';
            flow = 2; % 1 col, 2 row-wise normalization
            fvec=opt.Freqband(1):opt.Freqband(2);
            load(fullfile(Projfolder,animals{subj}),'srate');
            tvec = (tsec>=opt.TimeWin(1)) & (tsec<=opt.TimeWin(2));
            tsec = tsec(tvec);
            labels = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
            for roi = 1:2
                % (3) estimate MVAR coeffs using stok
                KF = dynet_SSM_STOK(epochs((1:6)+(roi-1)*6,(1:6)+(roi-1)*6,tvec),opt.ModOrds,ff);  
                % (4) estimate PDC based on MVAR coeffs
                PDC.(animID{subj})((1:6)+(roi-1)*6,(1:6)+(roi-1)*6,:,:) = dynet_ar2pdc(KF,srate,fvec,measure,keepdiag,flow);
            end
        end

        save(fullfile(Projfolder,SaveFileName),'PDC', 'fvec','tsec','labels');
    else
        load(fullfile(Projfolder,SaveFileName));
    end
    
    %% Normalize PDC values using prestim PDCs
    if opt.NormPrestim
        IndNorm = find(tsec<0);
        IndNorm = IndNorm(100:end); % to remove the unstable part of the PDCs
         for subj = 1:numel(animals)
            PDC.(animID{subj}) = (PDC.(animID{subj})-mean(PDC.(animID{subj})(:,:,:,IndNorm),4));
            
         end
         
         SaveFigName = [SaveFileName '_NormPrestim'];
    else
        SaveFigName = SaveFileName;
    end

    %% Plot the results: individuals and average layer connectivities
    ROIs = {'cS1','iS1'};
    
    if opt.plotfig
        if ~exist(opt.figpath,'dir')
            mkdir(opt.figpath);
        end
        
        for subj = 1:numel(animals)
            Data = PDC.(animID{subj})(:,:,:,:);
            % what kind of normalization to be done
            if strcmpi(opt.Normalize,'Channel')% (1) normalize over every single channel
                Data = Data./repmat(max(max(Data,[],4),[],3),[1 1 size(Data,3) size(Data,4)]); 
                
            elseif strcmpi(opt.Normalize,'All')% (2) normalize over all channels
                for ch = 1:size(Data,1)
                    Data(ch,ch,:,:) = 0;
                end
                Data = Data./max(Data(:));
            end
            % base-line correction? Does it make sense here?
            Databm  =   repmat(mean(Data(:,:,:,:),4),[1 1 1 size(Data,4)]);
            Databs  =   repmat(std(Data(:,:,:,:),[],4),[1 1 1 size(Data,4)]);
            Datap   =   Data; %(Data-Databm)./Databs;
            %plot individual animal resaults
            if false
                FIG = dynet_connplot(Datap(1:6,1:6,:,:),tsec,fvec,labels, [], [], [],0);
                set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
                export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_' animID{subj} '_cS1']),'-pdf');

                FIG = dynet_connplot(Datap(7:12,7:12,:,:),tsec,fvec,labels, [], [], [],0);
                set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
                export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_' animID{subj} '_cS1']),'-pdf');
            end
            close all;
            
            DataM(:,:,:,:,subj) = Data;%(Data-Databm)./Databs;
        end
        
        % plot average over all animals
        if opt.NormPrestim
            ranges = [-1 1];
        else
            ranges = [0 1];
        end
        TimeInd = (tsec>=-50);
        Data = (mean(DataM(1:6,1:6,:,TimeInd,:),5));
        FIG = dynet_connplot(Data./max(Data(:)),tsec(TimeInd),fvec,labels, ranges, [], [],0);
        if opt.NormPrestim, colormap(jmaColors('coolhotcortex'));end
        set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
        export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_ AverageAll_' ROIs{1}]),'-pdf');
        
        Data = mean(DataM(7:12,7:12,:,TimeInd,:),5);
        FIG = dynet_connplot(Data./max(Data(:)),tsec(TimeInd),fvec,labels, ranges, [], [],0);
        if opt.NormPrestim, colormap(jmaColors('coolhotcortex'));end
        set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
        export_fig(FIG,fullfile(opt.figpath,[SaveFigName '_ AverageAll_' ROIs{2}]),'-pdf');
        
              
    end
    
    %% estimate UPWARD and DOWNWARD FCs
    clear DataM
    for subj = 1:numel(animals)% Aggregate animals data
        
        Data = PDC.(animID{subj})(:,:,:,:);
        for ch = 1:size(Data,1)
            Data(ch,ch,:,:) = 0;
        end
        Data = Data./max(Data(:));      
        DataM(:,:,:,:,subj) = Data;
    end 
    
    Data = mean(DataM(:,:,:,:,:),5);
   
    %% (1) node-wise analysis
    load('LayerColors.mat');
    LNames = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
    %-----------------------Nodes outflows-----------------------------------
    if true
        lnum  = size(Data,1)/2;
        for roi = 1:2
            Fig1 = figure;
            for i = 1:lnum
                subplot(lnum,1,i);
                for j  = 1:lnum
                    SP(j)=plot(tsec,squeeze(mean(Data(j+(roi-1)*lnum,i+(roi-1)*lnum,:,:),3)),'color',Colors(j,:),'linewidth',2);
                    hold on;
                end
                title(['L' num2str(i)]);
                xlim([-50 100])
                if opt.NormPrestim
                    ylim([-0.4 0.5])
                else
                    ylim([0 .7])
                end
                vline([0],'k--')
            end
            legend(SP,LNames);
            xlabel('Time (ms)');
            ylabel('outPDC')
            set(Fig1,'unit','inch','position',[2 0 15 20],'color','w')
            export_fig(Fig1,fullfile(opt.figpath,[SaveFigName '_LayersOutflow_' ROIs{roi}]),'-pdf');
        end
    end

    %% -----------------------Average Outflows---------------------------------
    if true
        for roi = 1:2
            lnum  = size(Data,1)/2;
            Data2 = Data((1:lnum)+(roi-1)*lnum,(1:lnum)+(roi-1)*lnum,:,:);
            PDCavgOut = arrayfun(@(x) squeeze(mean(Data2([1:x-1 x+1:end],x,:,:),1)),1:size(Data2,1),'uni',false);

            Fig2 = figure;
            line([-100 300],[0 0],'linestyle','--','color','k','linewidth',1.3);
            hold on;
            for i = 1:numel(PDCavgOut)
                SP(i) = plot(tsec,mean(PDCavgOut{i},1),'color',Colors(i,:),'linewidth',2);
            end
            xlim([-50 100])
            xlabel('Time (ms)')
            vline([0],'k--')
            legend(SP,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false));
            set(gca,'fontsize',16);
            set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
            title ('Average Outflow of layers')
            export_fig(Fig2,fullfile(opt.figpath,[SaveFigName '_AveragewOutflow_' ROIs{roi}]),'-pdf');
        end
    end
    %close all;
    %% ---------------------------Directionality-----------------------------
    FS = 12;
    
    for roi =1:2
       lnum  = size(DataM,1)/2;
       DataM2 = abs(DataM((1:lnum)+(roi-1)*lnum,(1:lnum)+(roi-1)*lnum,:,:,:)); 
       %Diff = arrayfun(@(x) squeeze(mean(DataM2(1:x-1,x,:,:,:),1) - mean(DataM2(x+1:end,x,:,:,:),1))./squeeze(mean(cat(1,mean(DataM2(1:x-1,x,:,:,:),1),mean(DataM2(x+1:end,x,:,:,:),1)),1)),1:6,'uni',false);
       Diff = arrayfun(@(x) squeeze(mean(DataM2(1:x-1,x,:,:,:),1) - mean(DataM2(x+1:end,x,:,:,:),1)),1:6,'uni',false);%./squeeze(mean(Data(x,x:end,:,:),2) + mean(Data(x,1:x,:,:),2)),1:size(Data,1),'uni',false);
       Diff = cellfun(@(x) mean(x,3),Diff,'uni',false);
       
       if true
           %-------------------------Averaged Over Frequencies-----------------
            Fig2 = figure;
            line([-50 100],[0 0],'linestyle','--','color','k','linewidth',1.3);
            hold on;
            for i = 2:numel(PDCavgOut)-1
                SP(i-1) = plot(tsec,sum(Diff{i},1),'color',Colors(i,:),'linewidth',2);
            end
            xlim([-50 100])
            xlabel('Time (ms)')
            vline([0],'k--')
            legend(SP,arrayfun(@(x) ['L' num2str(x)],2:5,'uni',false));
            set(gca,'fontsize',16);
            set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
            title ('Upward flow of layers')
            export_fig(Fig2,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AveragedFreqs_' ROIs{roi}]),'-pdf');

            %--------------------------All Frequencies-------------------------
            ind = (tsec>=-100);
            tseccorr = tsec(ind);
            Fig22 = figure;
            for i = 2:numel(PDCavgOut)-1
                SP(i-1) = subplot(4,1,i-1);
                Dataplot = Diff{i}(:,ind);
                imagesc(Dataplot);
                axis xy;
                xtickso = find(ismember(tseccorr,-100:100:300));
                set(gca,'xtick',xtickso,'xticklabel',tseccorr(xtickso),'ytick',1:40:numel(fvec),'yticklabel',fvec(1:40:end));
                M = max(abs(Diff{i}(:)));
                %caxis([-M M]);
                caxis([-.2 .2]);
                vline(xtickso(2),'k--');
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
            DiffM = mean(cat(3,Diff{2:5}),3);
            Fig3 = figure;
            imagesc(DiffM(:,ind));
            axis xy;
            xtickso = find(ismember(tseccorr,-100:100:300));
            set(gca,'xtick',xtickso,'xticklabel',tseccorr(xtickso),'ytick',1:40:numel(fvec),'yticklabel',fvec(1:40:end));
            caxis([-.07 .07]);
            vline(xtickso(2),'k--');
            colorbar;
            xlabel('Time(msec)');
            ylabel('Frequency(Hz)');
            set(gca,'fontsize',FS);
            colormap(jmaColors('coolhotcortex'))
            set(Fig3,'unit','inch','position',[1 1 10 4],'color','w');
            export_fig(Fig3,fullfile(opt.figpath,[SaveFigName '_Upwardflow_Norm_AveragedLayers_' ROIs{roi}]),'-pdf');

            %close all
       end

    end
    
 %% =================================------Second part: Statistics-----===============================  
 
 
%%    Next step is to test the significance of up/down directionalities: Using permutation test 
if opt.doStats      
        LNum = size(Data,1);
        if strcmpi(opt.typePermute,'layers')% Permute layers
            Perms = perms(1:LNum/2);% Generate all possible permutations
            Nperm = size(Perms,1);
        elseif strcmpi(opt.typePermute,'all')% Permute all FC values
            Nperm = opt.NPerm;
        end
        
        PermFileName = fullfile(Projfolder,['PDCSTOK_PermuteState_' opt.typePermute '_' opt.SStat '.mat']);
        if ~exist(PermFileName,'file') || opt.doPermute
            %tic
            for p = 1:Nperm % for each permutation:
                if mod(p,3)==0, disp(p); end
                if strcmpi(opt.typePermute,'layers')% Permute layers
                    Ps = Perms(p,:);
                    Ps = [Ps Ps+LNum/2];
                    Datap = Data(Ps,Ps,:,:); 
                elseif strcmpi(opt.typePermute,'all')% Permute all FC values
                    if p~=Nperm
                        Datap = PermuteElements(Data,LNum);
                    else
                        Datap = Data;
                    end
                end
                % Calculate directionality of FCs, and generate the null (random) distribution
                for t = 16:size(Data,4)
                    for f = 1:size(Data,3)
                        for l =1:2
                            Up  = triu(squeeze(Datap((1:6)+(LNum/2*(l-1)),(1:6)+(LNum/2*(l-1)),f,t)),1);
                            Down  = tril(squeeze(Datap((1:6)+(LNum/2*(l-1)),(1:6)+(LNum/2*(l-1)),f,t)),-1);
                            if strcmpi(opt.SStat,'Mean')
                                UpwardN(f,t,p,l) = (sum(Up(:)) - sum(Down(:)))/(sum(Up(:)) + sum(Down(:))); % estimate the mean difference
                            elseif strcmpi(opt.SStat,'Paired')
                                Up = Up';
                                [~,~,ci,stats]= ttest(Up(Up~=0),Down(Down~=0),'Tail','both');
                                UpwardN(f,t,p,l) = stats.tstat;
                                DF = stats.df;
                            end
                        end
                    end
                end
                %toc
            end

            % Indicate p-values for each point in time and frequency, based
            % on the null distribution
            for p = 1:Nperm
                if mod(p,20)==0, disp(p); end
                for t = 1:size(Data,4)
                    for f = 1:size(Data,3)
                        for l = 1:2
                            PvalsU(f,t,p,l) = sum(UpwardN(f,t,:,l)>=UpwardN(f,t,p,l))/Nperm;
                            PvalsD(f,t,p,l) = sum(UpwardN(f,t,:,l)<=UpwardN(f,t,p,l))/Nperm;
                        end
                    end
                end
            end
            save(PermFileName,'UpwardN','PvalsU','PvalsD');
        else
            load(PermFileName,'UpwardN','PvalsU','PvalsD');
        end 
        
        % extract cluster statistics
%         Suprath = .05;% indicate suprathreshold
%         for p = 1:Nperm-1
%             CS= findclust(PvalsU(:,:,p),UpwardN(:,:,p),Suprath);
%             ClusterStats(p) = CS(1);
%         end
%         [CS Clusters]= findclust(PvalsU(:,:,end),UpwardN(:,:,end),Suprath);
    
    
    %% visualize results
     % (2) general
    
%     for t = 1:size(Data,4)
%         for f = 1:size(Data,3)
%             Up  = triu(squeeze(Data(:,:,f,t)),1);
%             Down  = tril(squeeze(Data(:,:,f,t)),-1);
%             Upward(f,t) = sum(Up(:)) - sum(Down(:));
%             UpwardN(f,t) = (sum(Up(:)) - sum(Down(:)))/(sum(Up(:)) + sum(Down(:)));
%         end
%     end

    ROIs = {'cS1','iS1'};
    for l = 1:2
        FIG = figure;
        %UpwardN(abs(UpwardN)<.05)=0;
        imagesc(UpwardN(:,:,end,l),'alphadata',((PvalsU(:,:,end,l)<0.05)+(PvalsD(:,:,end,l)<0.05))*.8+.2); axis xy;
        %imagesc(UpwardN(:,:,end,l));
        % Set figure designations
        tTicks = -100:20:100;
        [~,Ind]=ismember(tTicks,tsec);
        vline(Ind(6),'k--')
        set(gca,'xtick',Ind,'xticklabel',tsec(Ind));
        axis xy
        %colormap(jmaColors('coolhotcortex'));
        colormap('jet');
        caxis([-3 3])    
        xlabel('Time(mSec)');
        ylabel('Frequency(Hz)');
        %xlim([20 751])

        colorbar;
        text(760,255,'Upward','color','r','fontweight','bold')
        text(760,-5,'Downward','color','b','fontweight','bold')
        set(FIG,'unit','inch','position',[1 4 15 5])
        export_fig(FIG,fullfile(opt.figpath,['UPDOWN_STOKPDC_AverageAll_permcorrection' opt.Normalize '_' ROIs{l} '_' opt.SStat]),'-pdf');
    end    

    
    %% What about bootstrapping
    
    % calculating bootstrapping in R!!!
    % For that pupose use R CMD Batch with 'system' command in matlab
    % check this out : https://www.youtube.com/watch?v=Guw2XgGvl44
    % https://stackoverflow.com/questions/6695105/call-r-scripts-in-matlab
    
end    
end

%% Extra functions
function Datap = PermuteElements(Data,LNum)
%  this function permutes Laminar LFP data, where Data is a
%  LNUMxLNUMxFreqxTime matrix
    RI = randperm(LNum*(LNum-1),LNum*(LNum-1));
    OI = 1:30;
    OIorig = reshape(1:LNum*LNum,LNum,LNum); 
    RIorig = OIorig;
    ind = 1;
    for i = 1:LNum
        for j = 1:LNum
            if i~=j
                IndR(j,i) = RI(ind);
                IndO(j,i) = OI(ind);
                ind = ind+1;       
            end
        end
    end
    for i = 1:numel(IndR)
        RIorig(find(IndO==i)) = OIorig(find(IndR==i));
    end

    Datat = reshape(Data,LNum*LNum,size(Data,3),size(Data,4));
    Datap = reshape(Datat(reshape(RIorig,LNum*LNum,1),:,:),LNum,LNum,size(Data,3),size(Data,4)); 
end

