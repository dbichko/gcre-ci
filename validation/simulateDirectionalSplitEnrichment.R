# sim data
require(tidyverse)

caseNo = 50
controlNo = 100
numGenes = 5

maxIter = 100
maxTrials = 10

caseIdx = 1
controlIdx = caseIdx + caseNo

genPatMat = function(freqOnes, flipVec) {
  patMat = matrix(0, nrow = 1 + numGenes, ncol = caseNo + controlNo)
  for (i in 1:numGenes) {
    if (flipVec[i]==1) {
      patMat[i+1,sample(1:caseNo, freqOnes*caseNo)] = 1
    } 
    if (flipVec[i]==-1) {
      patMat[i+1,sample(controlIdx:ncol(patMat), freqOnes*controlNo)] = 1
    }
    if (flipVec[i]==0) {
      patMat[i+1,sample(1:ncol(patMat), freqOnes*(caseNo+controlNo))] = 1
    }
  }
  patMat[1, ] = c(rep(0,caseNo),rep(1,controlNo))
  rownames(patMat) = c("pheno", 1:numGenes)
  patMat
}

enrichmentTest = function(patMat, flipVec, alternative="two") {
  if (is.vector(patMat)) return(-log10(0.5))
  genMat = patMat[-1,]
  if (nrow(patMat)>2) {
    andVec = colSums(genMat) > 0
  }
  else {
    andVec = genMat
  }
  -log10(fisher.test(andVec, patMat[1,], alternative=alternative)$p.value)
}

directionalSplitTestTwoSided = function(patMat, flipVec) {
  genMat = patMat[-1,]
  posMats = patMat[c(TRUE, flipVec>0),]
  negMats = patMat[c(TRUE, flipVec<0),]
  posP = enrichmentTest(posMats, flipVec)
  negP = enrichmentTest(negMats, flipVec)
  posP + negP
}

simulateIt = function(freqOnes, genFlipVec, testFlipVec, maxIter) {
  res = list()
  origPatMat = genPatMat(freqOnes, genFlipVec)
  for (i in 1:(maxIter+1)) {
    patMat = origPatMat
    if (i>1) {
      patMat = patMat[, sample(1:ncol(patMat))]
      patMat[1,] = origPatMat[1,]
    }
    enrich = enrichmentTest(patMat, testFlipVec)
    directionalSplit = directionalSplitTestTwoSided(patMat, testFlipVec)  
    res[[length(res)+1]] = data_frame(enrich=enrich, 
                                      directionalSplit = directionalSplit)
  }
  resAll = bind_rows(res)
  #print(resAll)
  pvals = (1 + sapply(resAll, 
                      function(col) { sum(col[-1] > col[1]) })) / maxIter
  data_frame(enrich=pvals[1], directionalSplit=pvals[2])
}

# error rate control with random data
res = list()
genFlipVec = c(0,0,0,0,0)
testFlipVec = c(1,1,1,1,1)
freqOnes = 0.5
for (trial in 1:maxTrials) {
  res[[length(res)+1]] = simulateIt(freqOnes, genFlipVec, testFlipVec, maxIter)
}
resAllRandom = bind_rows(res)


# all PLUS
res = list()
genFlipVec = c(1,1,1,1,1)
testFlipVec = c(1,1,1,1,1)
freqOnes = 0.5
for (trial in 1:maxTrials) {
  res[[length(res)+1]] = simulateIt(freqOnes, genFlipVec, testFlipVec, maxIter)
}
resAllPLUS = bind_rows(res)

# all MIXED
res = list()
genFlipVec = c(1,-1,1,-1,1)
testFlipVec = c(1,-1,1,-1,1)
freqOnes = 0.5
for (trial in 1:maxTrials) {
  res[[length(res)+1]] = simulateIt(freqOnes, genFlipVec, testFlipVec, maxIter)
}
resAllMIXED = bind_rows(res)


# wrong directions
res = list()
genFlipVec = c(1,-1,1,-1,1)
testFlipVec = c(-1,1,1,1,1)
freqOnes = 0.5
for (trial in 1:maxTrials) {
  res[[length(res)+1]] = simulateIt(freqOnes, genFlipVec, testFlipVec, maxIter)
}
resAllWRONG = bind_rows(res)


ggplot(resAllRandom, aes(x=enrich, y=directionalSplit)) + geom_point()
ggplot(resAllPLUS, aes(x=enrich, y=directionalSplit)) + geom_point()
ggplot(resAllMIXED, aes(x=enrich, y=directionalSplit)) + geom_point()
ggplot(resAllWRONG, aes(x=enrich, y=directionalSplit)) + geom_point()
