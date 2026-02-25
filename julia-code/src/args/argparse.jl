#= Simple argument parser for some CLI scripts =#
#/ Start module
module MArgParse

using ArgParse

#################
### FUNCTIONS ###
"Parse arguments for goodness-of-fit scripts"
function parsegof()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--numeps", "-n"
        help = "Number of ε to try when fitting heavy-tailed distribution(s)"
        arg_type = Int
        default = 100
        "--top", "-t"
        help = "Number of samples with the highest no. of reads to keep."
        arg_type = Int
        default = 10
        "--filter", "-f"
        help = "Flag for filtering data on, e.g., no. of samples, no. of components, etc."
        action = :store_true
    end
    return ArgParse.parse_args(s)
end

end # module MArgParse
#/ End module
