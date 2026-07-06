# Shared Taylor-law prediction-bin helpers for Figure 2B.
#
# The helper functions live in plot-fig2-A.jl. Keep this small wrapper so
# plot-fig2-B.jl can be included independently.

if !isdefined(Main, :Figure2A)
    Base.include(Main, joinpath(@__DIR__, "plot-fig2-A.jl"))
end

module TLPredictionBinPlotter

import Main.Figure2A:
    MM_TO_PT,
    NATURE_DOUBLE_WIDTH_PT,
    NATURE_AXIS_LABEL_PT,
    NATURE_TICK_PT,
    FONT_SCALE,
    LINE_WIDTH_SCALE,
    MARKER_SIZE,
    _normalize_prediction_columns!,
    _load_bin_data,
    _bin_point_groups,
    _add_icon!

end # module TLPredictionBinPlotter
