"""
    `streakpoke(r::DataFrames.DataFrameRow)`
    Given poke information as a DataFrame row compose a box whose y position is
    equal to the trial number, and its x extension is given by the poke duration
"""
function streakpoke(r::DataFrames.DataFrameRow;ordered = false)
    y_start = ordered ? r.order - 0.25 : r.Streak - 0.25
    y_finish = ordered ? r.order + 0.25 : r.Streak + 0.25
    x_start = r.In
    x_finish = r.Out
    Shape([(x_start,y_start),(x_start,y_finish),(x_finish,y_finish),(x_finish, y_start),(x_start,y_start)])
end

function outcome_col(r::DataFrames.DataFrameRow)
    if r.Reward
        fcol = :green
    elseif !r.Reward
        fcol = :grey
    end
    if r.Side != r.SideHigh
        fcol = :white
    end
    return fcol
end

function side_annotate!(p,r::DataFrames.DataFrameRow)
    if typeof(r.Side) <: Real
        side =  r.Side == 0 ? "R" : "L"
    else
        side = r.Side
    end
    annotate!(p,[(-0.5, r.Streak, Plots.text(side, 4, :center))])
end

"""
    `session_plot!(p,r::DataFrames.DataFrameRow)`
     Given poke information as a DataFrame row, plots a box on an existing
     figure with its positioned and color given by the trial number, poke time
     and outcome.
"""
function session_plot!(p,r::DataFrames.DataFrameRow;ordered = false)
    bp = streakpoke(r;ordered = ordered)
    fcol = outcome_col(r)
    plot!(p,bp, fillcolor = plot_color(fcol, 0.6),linecolor = :black,linewidth = 1)
    side_annotate!(p,r)
    p
end

"""
    `session_plot(df::DataFrames.AbstractDataFrame)`
     Given DataFrame containing all pokes in a session loops each row, to plot
     a box per each box.
"""

function session_plot(df::DataFrames.AbstractDataFrame; ordered = false)
    xyfont = font(18, "Bookman Light")
    legfont = font(14, "Bookman Light")
    xprop = ("Time (seconds)", (-1,16), xyfont)
    yprop = ("Trials", xyfont,false, false)
    plt = plot(legend = false, xaxis = xprop, yaxis = yprop, yticks = false)
    for r in eachrow(df)
        FLPDevelopment.session_plot!(plt, r;ordered = ordered)
    end
    plt
end
