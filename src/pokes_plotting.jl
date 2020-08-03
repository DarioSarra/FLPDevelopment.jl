function windowrange(dim, denom)
    range(0, 1, length = max(3, round(Int, dim/denom)))[2:end - 1]
end

function boxpoke(r::DataFrames.DataFrameRow)
    if typeof(r.Side) <: Real
        side =  r.Side
    else
        side = r.Side == "R" ? 0 : 1
    end
    y = (2 * side) - 1
    start = r.PokeIn
    finish = r.PokeOut
    Shape([(start,0),(start,y),(finish,y),(finish,0),(start,0)])
end

function outcome_col(r::DataFrames.DataFrameRow)
    if r.Reward
        fcol = :yellow
    elseif !r.Reward
        fcol = :grey
    end
    if r.Side != r.SideHigh
        fcol = :white
    end
    return fcol
end

function stim_col(r::DataFrames.DataFrameRow)
    if r.Stim == "true"
        lcol = :blue
    else
        lcol = :black
    end
    return lcol
end

function wall_scatter!(p,r::DataFrames.DataFrameRow, lcol)
    if  r.Wall == "true"
        y = (2 * r.Side) - 1
        start = r.PokeIn/1000
        finish = r.PokeOut/1000
        Width = finish - start
        I = windowrange(Width, 0.01)
        J = y.*[0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]
        pts = vec([(Float64(start + Width * i), Float64(1* j)) for i in I, j in J])
        scatter!(pts,
                markershape = :rect,
                markersize = 0.2,
                markeralpha = 0.6,
                markercolor = lcol,
                markerstrokealpha = 0,
                )
    end
end

function protocol_annotate!(p,r::DataFrames.DataFrameRow)
    if r.Delta!=0
        y = (((2 * r.Side) - 1) * 1.05)
        start = r.PokeIn/1000
        finish = r.PokeOut/1000
        Width = finish - start
        if r.Protocol == 0
            t = "L"
        elseif r.Protocol == 1
            t = "H"
        end
        x = start + (Width/2)
        annotate!([(x, y, Plots.text(t, 8, :center))])
    end
end

function poke_plot!(p,r::DataFrames.DataFrameRow)
    bp = boxpoke(r)
    lcol = stim_col(r)
    fcol = outcome_col(r)
    plot!(p,bp, fillcolor = plot_color(fcol, 0.6),linecolor = lcol,linewidth = 1)
    wall_scatter!(p,r,lcol)
    protocol_annotate!(p,r)
    p
end
