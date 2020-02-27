function reconMat = reduceMatrixRank(A,newRankVect)
[U, S, V] = svd(A);
reconMat = U(:,newRankVect)*S(newRankVect,newRankVect)*V(:,newRankVect)';
end