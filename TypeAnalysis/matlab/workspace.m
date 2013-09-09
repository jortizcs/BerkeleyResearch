plot(db1_7);
figure;
plot(db1_7(1:3600));
%% import data and take a slice
dataset = [temp_trace(1:3600), temp_trace(3601:7200), temp_trace(7201:10800)];
%dataset = [rh_trace(1:3600), rh_trace(3601:7200), rh_trace(7201:10800)];

%% try different k's for kmeans clustering
k=2;
[states, c, sumd, d] = kmeans(dataset(:,1), k);

gauss_params = [];
for s=1:k
    A1=dataset(find(states(:,1)==s), 1);
    gauss_params=[gauss_params; mean(A1) sqrt(var(A1))];
end

% ndist = makedist('Normal', 'mu', 1, 'sigma', 0.5);
% r = random(ndist);

%%  calculate the mean and variance for each class and instantiate
%   distribution to generate readings from accordingly

u_vals = unique(dataset);
A = [];
for i=1:3600
    A = [A, find(u_vals(:,1)==dataset(i,1))];
end
A = A';



%% learn the hmm and bin the unique values in the trace
[trans, emis] = hmmestimate(A, states);
[estTR,estE] = hmmtrain(A,trans,emis);

%%
[mseq, mstates] = hmmgenerate(3600, trans, emis);

mseq2=[];
for i=1:3600
    mseq2 = [mseq2; u_vals(mseq(1,i),1)];
end


%% test the classification accuracy of the classifier
[h, p, s]=kstest2(mseq2, rh_trace )

















