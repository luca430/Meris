#= Module to load arXiv dataset

=#
#/ Start module
module arXivLoader

#/ Packages
using Glob
using CSV, DataFrames, DataFramesMeta
using Random, StatsBase

#/ Modules, directories
import Meris.ARXIVDIR as ARXIVDIR

#/ STOPWORDS
const STOPWORDS = Set([
    "me", "my", "myself", "we", "our", "ours", "ourselves", "you", "your", "yours", "yourself", "yourselves", "he", "him", "his", "himself",
    "she", "her", "hers", "herself", "it", "its", "itself", "hey", "them", "their", "theirs", "themselves", "what", "which", "who", "whom",
    "this", "that", "these", "those", "am", "is", "are", "was", "were", "be", "been", "being", "have", "has", "had", "having", "do", "does", "did",
    "doing", "an", "the", "and", "but", "if", "or", "because", "as", "until", "while", "of", "at", "by", "for", "with", "about", "against", "between",
    "into", "through", "during", "before", "after", "above", "below", "to", "from", "up", "down", "in", "out", "on", "off", "over", "under", "again",
    "further", "then", "once", "here", "there", "when", "where", "why", "how", "all", "any", "both", "each", "few", "more", "most", "other", "some",
    "such", "no", "nor", "not", "only", "own", "same", "so", "than", "too", "very", "can", "will", "just", "don", "should", "now",
    "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"
])

#################
### FUNCTIONS ###
"Load all papers, put them into a single DataFrame"
function load(;
    DIR=ARXIVDIR * "processed/",
    stopwords=true
)

    # Make a dictionary with a structure that reflects the filesystem as in Luca's laptop:
    # domain -> topic -> samples (each sample = Vector{String} lines)
    domains = Dict{String,Dict{String,Vector{Vector{String}}}}()

    for domain_dir in filter(isdir, readdir(DIR; join=true))
        domain = splitpath(domain_dir)[end]
        topic_map = Dict{String,Vector{Vector{String}}}()

        for topic_dir in filter(isdir, readdir(domain_dir; join=true))
            topic = splitpath(topic_dir)[end]

            txt_files = filter(f -> endswith(f, ".txt"), readdir(topic_dir; join=true))
            if isempty(txt_files)
                continue
            end

            topic_map[topic] = [readlines(f) for f in txt_files]
        end

        if !isempty(topic_map)
            domains[domain] = topic_map
        end
    end

    # Create a DataFrame with words counts for all domains and topics
    big_df = DataFrame(
        domain=String[],
        topic=String[],
        component_id=String[],
        sample_id=String[],
        counts=Int[],
        nreads=Int[]
    )

    for (domain, topics) in domains
        for (topic, samples) in topics
            for (i, sample) in enumerate(samples)
                # Clean each word: keep only letters/numbers
                clean_sample = [replace(w, r"^[^A-Za-z0-9]+|[^A-Za-z0-9]+$" => "") for w in sample]

                # Remove empty or short words
                clean_sample = filter(!isempty, clean_sample)

                # Count words
                cnt_map = countmap(clean_sample)
                for (k, v) in cnt_map
                    push!(big_df, (domain, topic, k, domain[1:3] * topic * string(i), v, sum(values(cnt_map))))
                end
            end
        end
    end

    # Remove stopwords
    if !stopwords
        filter!(row -> !(row.component_id in STOPWORDS), big_df)
    end

    return big_df
end

end # module arXivSampler
#/ End module

