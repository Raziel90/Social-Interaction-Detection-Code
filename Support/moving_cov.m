function cov_seq=moving_cov(data,window)

if nargin<2
    window = 15;
end

% data = rand(378,9);
[M,N] = size(data);

numberCovarianceMatrices = M - window + 1;
covarianceMatrices = zeros(N,N,numberCovarianceMatrices);
for nc = 1:numberCovarianceMatrices
    covarianceMatrices(:,:,nc) = cov(data(nc:(nc+window-1),:));
end

cov_seq=covarianceMatrices;
end