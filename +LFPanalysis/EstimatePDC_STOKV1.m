function [PDC, fvec] = EstimatePDC_STOKV1(Projfolder,varargin)

%% set default values
    opt = ParseArgs(varargin, ...
        'IDs'       ,[],...
        'plotfig'  ,true, ...
        'ModOrds'  ,15,...
        'Cnd'      , 6, ...
        'recalc'   ,false, ...
        'Normalize','All',...
        'Freqband' , [5 150],...
        'doPermute', true,...
        'typePermute','all',...
        'TimeWin'  , [-200 300],...
        'figpath'  ,fullfile(fileparts(fileparts(Projfolder)),'Results') ...
        );

    animals = dir(fullfile(Projfolder,'*Layers.mat')); % files in the project folder
    animals = {animals.name};
    animID = cellfun(@(x) x(1:6),animals,'uni',false); % animal IDs
    
    if ~isempty(opt.IDs) % if a subset of animals are selected
        animIdx = strcmp(animID,opt.IDs);
        animID = animID(animIdx);
        animals = animals(animIdx);
    end
    
    opt.Freqband = round(opt.Freqband); % Round the frequency band for now
    
%% Estimate MVAR parameters using STOK algorithm
    xticklabels = -500:100:500;
    xticks = linspace (0,1250,numel(xticklabels));
    
    if ~exist(fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd) '.mat']),'file') || opt.recalc
        for subj = 1:numel(animals)
            
            %---------------------- Load required Data --------------------
            disp(animID{subj});
            LFP = load(fullfile(Projfolder,animals{subj}),'lfpMouse');
            load(fullfile(Projfolder,animals{subj}),'tsec');
            dims = size(LFP.lfpMouse);
            LFP.lfpMouse = reshape(LFP.lfpMouse, dims(1),dims(2),dims(3)*dims(4),dims(5)); %treat all exps (recordings) as one

            %--------------- Estimate MVAR params and PDC -----------------
            Cond = opt.Cnd;
            % (1) Prepare the data
            epochs = permute(LFP.lfpMouse(:,:,:,Cond), [3,2,1]); %trials, nodes, time
            % (2) prepare the parameters
            ff=.99;
            keepdiag = 1; % 
            measure='sPDC';
            %measure='PDCnn';
            flow = 2; % 1 col, 2 row-wise normalization
            fvec=1:150;
            load(fullfile(Projfolder,animals{subj}),'srate');
            tvec = (tsec>=opt.TimeWin(1)) & (tsec<=opt.TimeWin(2));
            tsec = tsec(tvec);
            labels = arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false);
            % (3) estimate MVAR coeffs using stok
            KF = dynet_SSM_STOK(epochs(:,:,tvec),opt.ModOrds,ff);  
            % (4) estimate PDC based on MVAR coeffs
            PDC.(animID{subj}) = dynet_ar2pdc(KF,srate,fvec,measure,keepdiag,flow);
        end

        save(fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd)]),'PDC', 'fvec','tsec','labels');
    else
        load(fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd)]));
    end

    %% Plot the results: individuals and average layer connectivities
    
    if opt.plotfig
        if ~exist(opt.figpath,'dir')
            mkdir(opt.figpath);
        end
        tvec = (tsec>=-100) & (tsec<=opt.TimeWin(2));
        tsec = tsec(tvec);
            
        for subj = 1:numel(animals)
            Data = PDC.(animID{subj})(:,:,opt.Freqband(1):opt.Freqband(2),tvec);
            % what kind of normalization to be done
            if strcmpi(opt.Normalize,'Channel')% (1) normalize over every single channel
                Data = Data./repmat(max(max(Data,[],4),[],3),[1 1 size(Data,3) size(Data,4)]); 
                
            elseif strcmpi(opt.Normalize,'All')% (2) normalize over all channels
                for ch = 1:size(Data,1)
                    Data(ch,ch,:,:) = 0;
                end
                Data = Data./max(Data(:));
            end
            % base-line correction? does it make sense here?
            % Databm = repmat(mean(Data(:,:,:,100:376),4),[1 1 1 size(Data,4)]);
            % Databs = repmat(std(Data(:,:,:,100:376),[],4),[1 1 1 size(Data,4)]);
            Datap = Data; %(Data-Databm)./Databs;
            %plot individual animal resaults
            FIG = dynet_connplot(Datap,tsec,fvec(opt.Freqband(1):opt.Freqband(2)),labels, [], [], [],0);
            set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
            %export_fig(FIG,fullfile(opt.figpath,['STOKPDC_' animID{subj} '_normalize_' opt.Normalize '_modelOrd_' num2str(opt.ModOrds)]),'-pdf');
            close all;
            
            DataM(:,:,:,:,subj) = Data;%(Data-Databm)./Databs;
        end
        
        % plot average over all animals
        Data = mean(DataM(:,:,:,:,:),5);
        FIG = dynet_connplot((Data./max(Data(:))),linspace(-100,opt.TimeWin(2),size(PDC.mID_40,4)),fvec(opt.Freqband(1):opt.Freqband(2)),labels, [], [], [],0);
        set(FIG,'unit','inch','position',[0 0 35 20],'color','w')
        export_fig(FIG,fullfile(opt.figpath,['STOKPDC_AverageAll_normalize_' opt.Normalize]),'-pdf');
        
        % generate a movie for average FC? (Maybe gamma [50 100])
        
    end
    
    %% estimate UPWARD and DOWNWARD FCs
    % Average the data
    clear DataM
    for subj = 1:numel(animals)-1 % On each animal separately
        
        Data = PDC.(animID{subj})(:,:,opt.Freqband(1):opt.Freqband(2),:);
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
    %-----------------------Nodes Inflows-----------------------------------
    if false
        Fig1 = figure;
        for i = 1:size(Data,1)
            subplot(size(Data,1),1,i);
            for j  = 1:size(Data,2)
                SP(j)=plot(tsec,squeeze(mean(Data(i,j,:,:),3)),'color',Colors(j,:),'linewidth',2);
                hold on;
            end
            title(['L' num2str(i)]);
            xlim([-100 300])
            ylim([0 0.16])
            vline2([0 37.5],{'k--','b--'})
        end
        legend(SP,LNames);
        xlabel('Time (ms)');
        set(Fig1,'unit','inch','position',[2 0 20 20],'color','w')
        export_fig(Fig1,fullfile(opt.figpath,['STOKPDC_NodesInflow' opt.Normalize]),'-pdf');
    end

    %% -----------------------Average Outflows---------------------------------
    PDCavgOut = arrayfun(@(x) squeeze(mean(Data([1:x-1 x+1:end],x,:,:),1)),1:size(Data,1),'uni',false);
    if false
        Fig2 = figure;
        line([-100 300],[0 0],'linestyle','--','color','k','linewidth',1.3);
        hold on;
        for i = 1:numel(PDCavgOut)
            SP(i) = plot(tsec,sum(PDCavgOut{i},1),'color',Colors(i,:),'linewidth',2);
        end
        xlim([-100 300])
        xlabel('Time (ms)')
        vline([0],'k--')
        legend(SP,arrayfun(@(x) ['L' num2str(x)],1:6,'uni',false));
        set(gca,'fontsize',16);
        set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
        title ('Average Outflow of layers')
        export_fig(Fig2,fullfile(opt.figpath,['STOKPDC_AveragewOutflow' opt.Normalize]),'-pdf');
    end
        %% -------------------------------Movie----------------------------------
    if true
        v = VideoWriter(fullfile(opt.figpath,['STOKPDC_MovieNet.avi']));
        v.FrameRate = 10;
        open(v);
        
        Nodesize = smoothdata(squeeze(mean(cat(3,PDCavgOut{:}),1))*80000);
        Stren = smoothdata(squeeze(mean(Data(:,:,:,:),3)),3,'movmean',10);
        Stren = Stren./max(Stren(:));
        
        FigM = figure;
        set(FigM,'unit','inch','position',[1 2 7 9],'color','w');
        Pos = [1.5 6;2.5 5; 2 4;1.5 2.; 2.5 1.5;2.1 .5];
         hold on;
        for T = 126:numel(tsec)

            for i = 1:6
                for j = 1:6
                    if i~=j
                        
                        if Stren(i,j,T)>0.25
                            if i>j
                                off = .1;
                            else
                                off = -.1;
                            end
                            if j==6
                                Li(i,j) = line([Pos(j,1)-.03 Pos(i,1)-.03],[Pos(j,2) Pos(i,2)]+off,'linewidth',Stren(i,j,T)*20,'color',[Colors(j,:) ]);
                            elseif i==6 && j==3
                                Li(i,j) = line([Pos(j,1)+.03 Pos(i,1)+.03],[Pos(j,2) Pos(i,2)]+off,'linewidth',Stren(i,j,T)*20,'color',[Colors(j,:) ]);
                            else
                                Li(i,j) = line([Pos(j,1) Pos(i,1)],[Pos(j,2) Pos(i,2)]+off,'linewidth',Stren(i,j,T)*20,'color',[Colors(j,:) ]);
                            end
                        end
                        %hold on;
                    end
                end
            end
            UP = sum(sum(triu(squeeze(mean(Data(:,:,:,T),3)))));
            Down = sum(sum(tril(squeeze(mean(Data(:,:,:,T),3)))));
            %A(1) = annotation('arrow',[.95 .95],[.2 .8],'linewidth',UP*20,'color','r') ;        
            %A(2) = annotation('arrow',[.05 .05],[.8 .2],'linewidth',Down*20,'color','b') ;

            N = scatter(Pos(:,1), Pos(:,2), Nodesize(T,:), Colors,'filled'); % nodes
            for l = 1:6, NT(l) = text(Pos(l,1)-.02,Pos(l,2),LNames{l},'Fontsize',14);end % node texts

            TT = text(1.8,6.5,['Time = ' num2str(tsec(T)) '(ms)'],'fontsize',14);

            xlim([1.2 2.8]);
            ylim([0 7])
            axis off;
            F(T-125) = getframe(FigM);
            pause(.01);
            delete(TT);
            delete(Li);
            delete(NT);
            delete(N);
            %delete(A)
        end
        
        %Mov = Movie(F,40);
        
        writeVideo(v,F);
        close(v);
        close all;
    end
    %% ---------------------------Directionality-----------------------------
    Diff = arrayfun(@(x) squeeze(mean(Data(1:x-1,x,:,:),1) - mean(Data(x+1:end,x,:,:),1))./squeeze(mean(cat(1,mean(Data(1:x-1,x,:,:),1),mean(Data(x+1:end,x,:,:),1)),1)),1:6,'uni',false);%./squeeze(mean(Data(x,x:end,:,:),2) + mean(Data(x,1:x,:,:),2)),1:size(Data,1),'uni',false);
    %Diff = arrayfun(@(x) squeeze(mean(Data(1:x,x,:,:),1) - mean(Data(x:end,x,:,:),1)),1:size(Data,1),'uni',false);
   if false
        Fig2 = figure;
        line([-100 300],[0 0],'linestyle','--','color','k','linewidth',1.3);
        hold on;
        for i = 2:numel(PDCavgOut)-1
            SP(i-1) = plot(tsec,sum(Diff{i},1),'color',Colors(i,:),'linewidth',2);
        end
        xlim([-100 300])
        xlabel('Time (ms)')
        vline([0],'k--')
        legend(SP,arrayfun(@(x) ['L' num2str(x)],2:5,'uni',false));
        set(gca,'fontsize',16);
        set(gcf,'unit','inch','position',[2 5 15 5],'color','w');
        title ('Upward flow of layers')
        export_fig(Fig2,fullfile(opt.figpath,['STOKPDC_Upwardflow' opt.Normalize '_Norm']),'-pdf');
   end
    
 %% =================================------Second part: Statistics-----===============================  
 
 
 
 %%    Next step is to test the significance of up/down directionalities: Using permutation test 
    
        LNum = size(Data,1);
        Perms = perms(1:LNum);% Generate all possible permutations
        Nperm = 2000;%size(Perms,1);
        PermFileName = fullfile(Projfolder,['PDCSTOK_Cnd' num2str(opt.Cnd) '_PermuteState_' opt.typePermute '.mat']);
        if ~exist(PermFileName,'file') || opt.doPermute
            for p = 1:Nperm % for each permutation:
                if mod(p,20)==0, disp(p); end
                
                if strcmp(opt.typePermute,'layers')% Permute layers
                    Datap = Data(Perms(p,:),Perms(p,:),:,:); 
                elseif strcmp(opt.typePermute,'all')% Permute all FC values
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
                % Calculate directionality of FCs, and generate the null (random) distribution
                for t = 1:size(Data,4)
                    for f = 1:size(Data,3)
                        Up  = triu(squeeze(Datap(:,:,f,t)),1);
                        Down  = tril(squeeze(Datap(:,:,f,t)),-1);
                        UpwardN(f,t,p) = (sum(Up(:)) - sum(Down(:)))/(sum(Up(:)) + sum(Down(:))); % estimate the mean difference
                    end
                end
            end

            % Indicate p-values for each point in time and frequency, based
            % on the null distribution
            for p = 1:Nperm
                if mod(p,20)==0, disp(p); end
                for t = 1:size(Data,4)
                    for f = 1:size(Data,3)
                        PvalsU(f,t,p) = sum(UpwardN(f,t,:)>=UpwardN(f,t,p))/Nperm;
                        PvalsD(f,t,p) = sum(UpwardN(f,t,:)<=UpwardN(f,t,p))/Nperm;
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
    FIG = figure;
    %UpwardN(abs(UpwardN)<.05)=0;
    imagesc(UpwardN(:,:,end),'alphadata',((PvalsU(:,:,end)<0.05)+(PvalsD(:,:,end)<0.05))*.8+.2); axis xy;
    % SET figure designations
    tTicks = -300:100:300;
    [~,Ind]=ismember(tTicks,tsec);
    vline(Ind(4),'k--')
    set(gca,'xtick',Ind,'xticklabel',tsec(Ind));
    axis xy
    %colormap(jmaColors('coolhotcortex'));
    colormap('jet');
    caxis([-0.8 .8]/2)    
    xlabel('Time(mSec)');
    ylabel('Frequency(Hz)');
    xlim([20 751])
    
    colorbar;
    text(760,155,'Upward','color','r','fontweight','bold')
    text(760,-5,'Downward','color','b','fontweight','bold')
    set(FIG,'unit','inch','position',[1 4 15 5])
    export_fig(FIG,fullfile(opt.figpath,['UPDOWN_STOKPDC_AverageAll_permcorrection' opt.Normalize]),'-pdf');
    
    
    %% What about bootstrapping
    
    % calculating bootstrapping in R!!!
    % For that pupose use R CMD Batch with 'system' command in matlab
    % check this out : https://www.youtube.com/watch?v=Guw2XgGvl44
    % https://stackoverflow.com/questions/6695105/call-r-scripts-in-matlab
    
    
end



