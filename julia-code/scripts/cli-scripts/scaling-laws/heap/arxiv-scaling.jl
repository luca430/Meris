#= Simple script to verify expected scaling of vocabulary with document size =#
#/ Packages
using DataFrames
using DataFramesMeta
using JLD2
using Random
using StatsBase

#/ Modules
using Meris
DATADIR = Meris.DATADIR * "heap/arxiv/"
mkpath(DATADIR)

#~ Specify variables
CATEGORY = "q-bio.PE"
save = true
npermutations = 64
rng  = Random.Xoshiro(npermutations)

#/ Get papers
#  note: currently gets papers from a specific category
arxivdf = Meris.arXivSampler.load_papers(; CATEGORY=CATEGORY)
#~ Check if requested number of permutations is possible
nsamples = length(arxivdf[!,:sample_id])
if nsamples < 9 && factorial(nsamples) > npermutations
    npermutations = factorial(nsamples)
end
#/ For each document permutation, compute the documentsize and vocabularysize
heapdf  = DataFrame(permutation=Int[], documentsize=Int[], vocabularysize=Int[])
for n in 1:npermutations
    vocabdf = Meris.arXivSampler.computevocabsize(arxivdf; rng=rng)
    vocabdf[!,:permutation] = fill(n, nrow(vocabdf))
    append!(heapdf, vocabdf)
end
#/ Compute also average scaling
#~ first, add a "row index" for each permutation
permutations = length(unique(heapdf[!,:permutation]))
heapdf[!,:rowidx] = repeat(1:nrow(heapdf) ÷ permutations, outer=permutations)
averagedf = @chain heapdf begin
    @by(:rowidx, :documentsize = mean(:documentsize), :vocabularysize=mean(:vocabularysize))
    @select(Not(:rowidx))
end

#~ Save
if save
    filename = "arxiv-vocabsize.jld2"
    jldsave(DATADIR*filename; raw = heapdf, average = averagedf)
end
