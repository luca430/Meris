#= Module to load arXiv dataset

=#
#/ Start module
module arXivLoader

#/ Packages
using CSV: makeunique
using Glob
using CSV, DataFrames, DataFramesMeta
using Random, StatsBase

#/ Modules, directories
import ..DataTools: filterdata
import Meris.ARXIVDIR as ARXIVDIR

#/ STOPWORDS
const STOPWORDS = Set([
    "me", "my", "myself", "we", "our", "ours", "ourselves",
    "you", "your", "yours", "yourself", "yourselves", "he", "him", "his", "himself",
    "she", "her", "hers", "herself", "it", "its", "itself",
    "hey", "them", "their", "theirs", "themselves", "what", "which", "who", "whom",
    "this", "that", "these", "those",
    "am", "is", "are", "was", "were", "be", "been", "being",
    "have", "has", "had", "having", "do", "does", "did",
    "doing", "an", "the", "and", "but", "if", "or", "because", "as",
    "until", "while", "of", "at", "by", "for", "with", "about", "against", "between",
    "into", "through", "during", "before", "after", "above", "below", "to", "from",
    "up", "down", "in", "out", "on", "off", "over", "under", "again",
    "further", "then", "once", "here", "there", "when", "where",
    "why", "how", "all", "any", "both", "each", "few", "more", "most", "other", "some",
    "such", "no", "nor", "not", "only", "own", "same", "so", "than", "too", "very",
    "can", "will", "just", "don", "should", "now",
    "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
    "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"
])

#################
### FUNCTIONS ###
"
    load()

Load *all* articles and put them into a single DataFrame. The DataFrame will be constructed so
that it can be used further in the analysis pipelines. It should contains [at least] the columns
    `component_id`, `sample_id`, `counts`, `nreads`
"
function load(
    ;
    DIR=ARXIVDIR * "processed/",
    minreads      = 4_000,
    mincomponents = 100_000,
    minsamples    = 30,
    stopwords     = true,
    applyfilter   = true,
    reorder       = true,
    top           = nothing,
    )
    #~ Allocate a dictionary as
    #  (arXiv domain) -> (topic) -> (samples [article])
    domaindict = Dict{String,Dict{String,Vector{Vector{String}}}}()

    #~ Loop through all subdirectories, load all articles, and store their information
    for DOMAINDIR in filter(isdir, readdir(DIR; join=true))
        domainname = splitpath(DOMAINDIR)[end]
        #~ Allocate another dictionary that collects the subdomain(s)
        subdomaindict = Dict{String,Vector{Vector{String}}}()

        for SUBDOMAINDIR in filter(isdir, readdir(DOMAINDIR; join=true))
            topic = splitpath(SUBDOMAINDIR)[end]
            
            TXTFILES = filter(f -> endswith(f, ".txt"), readdir(SUBDOMAINDIR; join=true))
            #~ Extract articles;
            #  here, ARTICLES contains a list of lists with the secondary list containing all
            #  the words of a single article, which can be filtered (when desired)
            ARTICLES = [readlines(f) for f in TXTFILES]
            
            #~ Filter articles that have insufficient total. no of words
            (applyfilter) && (filter!(article -> length(article) > minreads, ARTICLES))
            #     #~ Skip subdomains entirely if they do not sufficient distinct words
            #     totalcomponents = length(unique(reduce(vcat, ARTICLES)))
            #     (totalcomponents < mincomponents) && (continue)
            #     #~ Skip subdomains that have insuffient no. of articles after filtering
            #     (length(ARTICLES) < minsamples) && (continue)
            # end
            #~ Add the remaining articles to the subdomain dictionary
            subdomaindict[topic] = ARTICLES
        end
        #~ Add subdomain to domain dictionary
        if !isempty(subdomaindict)
            if applyfilter
                #~ Skip subdomains that have insuffient no. of articles after filtering
                articlelengths = map(x -> length(x), values(subdomaindict))
                (sum(articlelengths) < minsamples) && (continue)
                #~ Skip subdomains entirely if they do not sufficient distinct words
                bagsofwords = map(x -> reduce(vcat, x), values(subdomaindict))
                bagofwords = reduce(vcat, bagsofwords)
                totalcomponents = length(unique(bagofwords))
                (totalcomponents < mincomponents) && (continue)
            end
            domaindict[domainname] = subdomaindict
        end
    end

    # Create a DataFrame with words counts for all domains and topics
    df = DataFrame(
        class=String[],
        sample_id=String[],
        component_id=String[],
        counts=Int[],
        nreads=Int[]
    )

    #/ Loop through the domain dictionary and add components [words] for each topic
    for (domain, subdomains) in domaindict
        for (subdomain, articles) in subdomains
            for (i, article) in enumerate(articles)
                #~ Clean each word: keep only letters/numbers
                _article = [
                    replace(word, r"^[^A-Za-z0-9]+|[^A-Za-z0-9]+$" => "") for word in article
                        ]
                #~ Remove empty or short words
                _sample = filter(!isempty, _article)
                #~ Count words
                cm = countmap(_sample)
                #~ Push them into the DataFrame
                for (component_id, counts) in cm
                    push!(
                        df,
                        (
                            domain,
                            #~ Construct `sample_id` from domain, subdomain, and index
                            domain[1:3] * subdomain * string(i),
                            component_id,
                            counts,
                            sum(values(cm))
                        )
                    )
                end
            end
        end
    end

    if applyfilter
        #~ filter data
        df = filterdata(
            df; minsamples=minsamples, minreads=minreads, mincomponents=mincomponents,
            reorder=reorder, top=top
        )
    end

    #/ Finally, remove stopwords (if desired)
    (!stopwords) && (filter!(row -> !(row.component_id in STOPWORDS), df))
    return df
end

########################
### HELPER FUNCTIONS ###
"""
Print some stats of interest
These stats are mostly used within the context of manuscripts, presentations, etc.
"""
function printstats(df::DataFrame)
    ncomponents = nrow(df)
    narticles = length(unique(df.sample_id))
    ntopics = length(unique(df.topic))
    @info "Stats:" ncomponents narticles ntopics
end

end # module arXivLoader
#/ End module

