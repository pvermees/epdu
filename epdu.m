function epdu(data,errors,makeblocks,numsamples,ticksperbin,numiterations,coupled)
% epdu_old decouples numsamples from ticksperbin and numiterations
%
% an estimator of probability density and its uncertainties

load bimodal;
% bimodaldata = [10+rand(1,10)*10 70+rand(1,10)*5];
% bimodalerrors = 1+rand(1,20)*5;
data = bimodaldata;
errors = bimodalerrors;

coupled = 0;
makeblocks = 1; % if makeblocks == 0, plot as a smooth curve
                % if makeblocks == 1, plot as a histogram, 
                % if makeblocks == 2, plot both
% for each of these samples, generate ticksperbin x numiterations shifted histograms
ticksperbin = 1;
numiterations = 50;
numsamples = 50;
wannasmooth = 0; % get rid of sharp plateaus by gentle iterative smoothing
% coupled = 1; % == 1 if you want to systematically go through all steps of the average shifted histogram (ash)
padding = 2; % padding should be at least 1 

% generate numsamples random samples following the distribution of the measurement errors
random_data = generate_data_from_norm_data(data,errors,numsamples);
data_range = max(max(random_data)) - min(min(random_data));
numbins = calculatenumbins(data);
binwidth = data_range/numbins;
ageticks = linspace(min(min(random_data))-padding*binwidth,...
                    max(max(random_data))+padding*binwidth,...
                    (numbins+2*padding)*ticksperbin+1);

% initialize...
if coupled
    tickfactor = ticksperbin; % coupling the ash and the metropolis algorithm enlarges the model space by a tickfactor 
else
    tickfactor = 1;
end
hists = zeros(numsamples*tickfactor*numiterations,length(ageticks));
cdfs = zeros(numsamples*tickfactor*numiterations,length(ageticks));
for i=1:numsamples,
    [shs,scdfs] = sh(random_data(i,:),ageticks,ticksperbin,numiterations,coupled);
    hists((i-1)*tickfactor*numiterations+1:i*tickfactor*numiterations,:) = shs;
    cdfs((i-1)*tickfactor*numiterations+1:i*tickfactor*numiterations,:) = scdfs;
end

% calculate the alpha percentiles
alpha = 5;
% leave out the last column of the hists and cdfs because of the way histc (in function sh) works
[hist_05,hist_95,hist_50,cdf_05,cdf_95,cdf_50] = percentiles(hists,cdfs,alpha);

% plot the results
titel = ['                                               ' ...
         ' # samples = ' num2str(numsamples)  ...
         ', # ticks/bin = ' num2str(ticksperbin)  ...
         ', # iterations = ' num2str(numiterations)     ];
plotresults(data,errors,hist_05(1:end-1),hist_95(1:end-1),hist_50(1:end-1),...
            cdf_05(1:end-1),cdf_95(1:end-1),cdf_50(1:end-1),ageticks,makeblocks,titel,wannasmooth);
            % omit the last entry because it is zero anyway

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
function random_data = generate_data_from_norm_data(data,errors,numsamples)
% function generate_data_from_norm_data(data,errors,numsamples)
% generates numsamples samples from the normally distributed
% data with standard deviation = errors.
% data and errors must be row vectors
if (numsamples==1), % if we want just histograms, we should take the expected values of the samples
    random_data = data;
else
    randnums = randn(numsamples,length(data));
    for i=1:numsamples,
        random_data(i,:) = randnums(i,:).*errors + data;
    end
end

function numbins = calculatenumbins(data)
numbins = round(1+log2(length(data))); % for now, Sturge's rule for number of bins

function [shs,scdfs] = sh(data,dataticks,ticksperbin,numiterations,coupled)
% shifted histograms with uncertainties
binwidth = dataticks(2)-dataticks(1);
numticks = length(dataticks);
histbase = zeros(size(dataticks));
if coupled
    iterator = [1:ticksperbin]; % go through all the 
else
    iterator = ceil(rand*ticksperbin); % choose a random offset
end
for i=iterator,
    edgeindexes = [i:ticksperbin:numticks-ticksperbin+i];
    edges = dataticks(edgeindexes);
    histvalues = histc(data,edges);
    hists = confidence_hist(histvalues,numiterations); % this will generate numiterations histograms
    for j=1:numiterations, % add these histograms to the set
        histogram = histbase;
        cumdistfn = histbase;
        for k=1:length(edgeindexes)-1,
            histogram(edgeindexes(k):edgeindexes(k+1)-1) = hists(j,k);
        end
        cumdistfn = cumsum(histogram,2)/ticksperbin; % if you don't divide, you count the same bins ticksperbin times
        shs((i-1)*numiterations*coupled+j,:) = histogram/sum(histogram*binwidth,2); % normalize area under histogram
        scdfs((i-1)*numiterations*coupled+j,:) = cumdistfn;
    end
end

function models = confidence_hist(histvalues,numiterations)
% construct numiterations histograms from the histvalues.
numgrains = sum(histvalues);
numbins = length(histvalues);
models = histvalues; % initialize the histogram to the true one
likelihoods(1) = likelihood(models(1,:),histvalues,numgrains);
i=2;
while (i <= numiterations),
    % initialized new model:
    models(i,:) = models(i-1,:);
    % randomly choose 2 bins:
    modifiedbins = ceil(rand(1,2)*numbins);
    while modifiedbins(1)==modifiedbins(2),
        modifiedbins = ceil(rand(1,2)*numbins); % make sure you chose 2 different bins!
    end
    % remember the sum of these two bins:
    sumofbins = sum(models(i-1,modifiedbins));
    % randomly choose new values between 0 and sumofbins:
    firstbin = round(rand*sumofbins);
    secondbin = sumofbins - firstbin;
    % create the new model:
    models(i,modifiedbins) = [firstbin,secondbin];
    % calculate the likelihood:
    likelihoods(i) = likelihood(models(i,:),histvalues,numgrains);
    if likelihoods(i) >= likelihoods(i-1)
        i = i+1;
    else
        P = likelihoods(i)/likelihoods(i-1);
        p = rand;
        if (p < P)
            i = i+1;
        end
    end  
end

function L = likelihood(n_model,n_obs,N)
    L = prod([binopdf(n_model(n_obs>0 & n_obs<N),N,n_obs(n_obs>0 & n_obs<N)./N) ...
          (1-n_model(n_obs==0)./N).^N (n_model(n_obs==N)./N).^N]);

function [hist_l,hist_u,hist_50,cdf_l,cdf_u,cdf_50] = percentiles(hists,cdfs,alpha)
if size(hists,1)==1 % if there is only one row, calculating percentiles doesn't make sense
    hist_u = hists;
    hist_50 = hists;
    hist_l = hists;
    cdf_u = cdfs;
    cdf_50 = cdfs;
    cdf_l = cdfs;    
else
    % calculate the percentiles of the columns of the hists and cdfs
    hist_u = prctile(hists,100-alpha/2);
    hist_50 = prctile(hists,50);
    hist_l = prctile(hists,alpha/2);
    cdf_u = prctile(cdfs,100-alpha/2);
    cdf_50 = prctile(cdfs,50);
    cdf_l = prctile(cdfs,alpha/2);
end

function plotresults(data,errors,hist_l,hist_u,hist_50,cdf_l,cdf_u,cdf_50,dataticks,makeblocks,titel,wannasmooth)
% if makeblocks == 1, then you plot the data as histograms,
% if makeblocks == 0, then you plot the histograms as curves (centered at middle of bins)
figure;
subplot(1,2,1);
if (makeblocks==0 | makeblocks==2)
    % plot the kernel estimate
    [x,fhat] = prob_dens_plot(dataticks,data,errors,100);
    plot(x,fhat,'linewidth',4,'color',[0.75 0.75 0.75]); hold on;
end
% plot the epdu histogram
if (makeblocks==1 | makeblocks==2)
    [blockhist_u,blockticks] = blockify(hist_u,dataticks);
    [blockhist_50,blockticks] = blockify(hist_50,dataticks);
    [blockhist_l,blockticks] = blockify(hist_l,dataticks);
    plot(blockticks,blockhist_u,'-k'); hold on;
    plot(blockticks,blockhist_50,'-k','linewidth',2); hold on;
    plot(blockticks,blockhist_l,'-k'); hold on;
end
if (makeblocks == 0 | makeblocks == 2)
    smoothdataticks = dataticks(1:end-1) + diff(dataticks)/2; % take the midpoints
    if (wannasmooth)
        hist_u_smooth = smooth(smoothdataticks,hist_u);
        hist_50_smooth = smooth(smoothdataticks,hist_50);
        hist_l_smooth = smooth(smoothdataticks,hist_l);
    else
        hist_u_smooth = hist_u;
        hist_50_smooth = hist_50;
        hist_l_smooth = hist_l;
    end
    plot(smoothdataticks,hist_u_smooth,'-k'); hold on;
    plot(smoothdataticks,hist_50_smooth,'-k','linewidth',2); hold on;
    plot(smoothdataticks,hist_l_smooth,'-k'); hold on;
end
title(titel);
xlabel('X');
ylabel('probability density');
%xlim([max(0,dataticks(1)) dataticks(end)]);
subplot(1,2,2);
% plot kernel estimate
if (makeblocks==0 | makeblocks==2)
    plot(x,length(data)*cumsum(fhat)/sum(fhat),'linewidth',4,'color',[0.75 0.75 0.75]); hold on;
end
if (makeblocks==1 | makeblocks==2)
    [blockcdf_u,blockticks] = blockify(cdf_u,dataticks);
    [blockcdf_50,blockticks] = blockify(cdf_50,dataticks);
    [blockcdf_l,blockticks] = blockify(cdf_l,dataticks);
    plot(blockticks,blockcdf_u,'-k'); hold on;
    plot(blockticks,blockcdf_50,'-k','linewidth',2); hold on;
    plot(blockticks,blockcdf_l,'-k'); hold on;
end
if (makeblocks == 0 | makeblocks == 2)
    if (wannasmooth)
        cdf_u_smooth = smooth(smoothdataticks,cdf_u);
        cdf_50_smooth = smooth(smoothdataticks,cdf_50);
        cdf_l_smooth = smooth(smoothdataticks,cdf_l);
    else
        cdf_u_smooth = cdf_u;
        cdf_50_smooth = cdf_50;
        cdf_l_smooth = cdf_l;
    end
    plot(smoothdataticks,cdf_u_smooth,'-k'); hold on;
    plot(smoothdataticks,cdf_50_smooth,'-k','linewidth',2); hold on;
    plot(smoothdataticks,cdf_l_smooth,'-k'); hold on;
end
xlabel('X');
ylabel('cumulative # observations');
%xlim([max(0,dataticks(1)) dataticks(end)]);

function smooth_data = smooth(ticks,data)
% interpolate data that have plateaus in them
smooth_data = data;
[start,finish] = findplateauedges(data);
if (start~=0) % if there are plateaus
    for i=1:length(start),
        if (start(i)==1 & finish(i)==length(data)) % pathetic dataset
            return
        elseif (start(i)==1)
            plateauindexes = [start(i):finish(i)+1];
            smooth_data(plateauindexes) = leftsmooth(smooth_data(plateauindexes),ticks(plateauindexes));
        elseif (finish(i)==length(data))
            plateauindexes = [start(i)-1:finish(i)];
            smooth_data(plateauindexes) = rightsmooth(smooth_data(plateauindexes),ticks(plateauindexes));
        else
            plateauindexes = [start(i)-1:finish(i)+1];
            smooth_data(plateauindexes) = recursivesmooth(smooth_data(plateauindexes),ticks(plateauindexes));
        end
    end
end

function smooth_data = leftsmooth(data,ticks)
plateaulength = length(data)-2;
smooth_data = data;
if (plateaulength==0)
    return
end
smooth_data(end-1) = interp1([ticks(end-2) ticks(end)],[data(end-2) data(end)],ticks(end-1));
smooth_data(1:end-1) = leftsmooth(smooth_data(1:end-1),ticks(1:end-1));

function smooth_data = rightsmooth(data,ticks)
plateaulength = length(data)-2;
smooth_data = data;
if (plateaulength==0)
    return
end
smooth_data(2) = interp1([ticks(1) ticks(3)],[data(1) data(3)],ticks(2));
smooth_data(2:end) = rightsmooth(smooth_data(2:end),ticks(2:end));

function smooth_data = recursivesmooth(data,ticks)
plateaulength = length(data)-3;
smooth_data = data;
if (plateaulength==0)
    return
end
if (plateaulength==1)
    if (data(1)<data(2) & data(3)<data(4)) | ...
       (data(1)>data(2) & data(3)>data(4)) % monotone curve
        smooth_data(2:3) = interp1([ticks(1) ticks(4)],[data(1) data(4)],ticks(2:3));
    end
    % otherwise you're on a peak or in a valley, and you shouldn't change the data:
    return
end
smooth_data(2) = interp1([ticks(1) ticks(3)],[data(1) data(3)],ticks(2));
smooth_data(end-1) = interp1([ticks(end-2) ticks(end)],[data(end-2) data(end)],ticks(end-1));
smooth_data(2:end-1) = recursivesmooth(smooth_data(2:end-1),ticks(2:end-1));

function [plateaustartindexes,plateauendindexes] = findplateauedges(data)
% returns [0,0] if there are no plateaus
plateauindexes = find(diff(data)==0); % the locations of the plateaus
if (isempty(plateauindexes))
    plateaustartindexes = 0;
    plateauendindexes = 0;
    return
end
indexindexes = find(diff(plateauindexes)~=1)+1;
plateaustartindexes = plateauindexes([1 indexindexes]);
plateauendindexes = [plateauindexes(indexindexes-1)+1 plateauindexes(end)+1];

function [blockhist,blockticks] = blockify(histogram,blockedges)
% dataticks has to be a row vector
blockticks = sort([blockedges blockedges]);
blockhist(1) = histogram(1);
blockhist(length(blockticks)) = histogram(end);
blockhist(2:2:length(blockticks)-2) = histogram;
blockhist(3:2:length(blockticks)-1) = histogram;

function [x,fhat] = prob_dens_plot(x,X,sigma,resolution)
% function fhat = prob_dens_plot(X,sigma,resolution) based on Silverman (1987)
% evaluate the data X at the values of x
numsamples = length(X);
x = linspace(min(x),max(x),resolution);
fhat = zeros(size(x));
for i=1:length(x),
    fhat(i) = f(x(i),X,sigma,numsamples);
end

function fhat = f(x,X,sigma,numsamples)
% x should be a scalar, not a vector!
fhat = sum(normpdf(x,X,sigma))/numsamples;