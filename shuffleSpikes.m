function shuffledSpikes = shuffleSpikes(spTrain)
spPrime = diff(spTrain);
shuffleIdx = randperm(numel(spPrime));
shuffledSpikes = cumsum([spTrain(1),spPrime(shuffleIdx)]);
end