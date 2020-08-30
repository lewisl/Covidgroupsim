#########################################################################################
# tracking.jl
#########################################################################################

pyplot()   # initialize plotting backend for Plots


"""
- use incr!(ctr, :day) for day of the simulation:  creates and adds 1
- use reset!(ctr, :day) to remove :day and return its current value, set it to 0
- use ctr[:day] to return current value of day
"""
const ctr = counter(Symbol) # from package DataStructures

# for debugging simulations: daily outcome entries as named tuples
const spreadq = []
const transq = []
const tntq = []
const r0q = []


# for Johns Hopkins US actual data
struct Col_ref
    date::String
    col::Int64
end


# tracking statistics

function showq(qname)
    for item in qname
        println(item)
    end
end


function reviewdays(q=spreadq)
    for it in q
        println(it)
        print("\nPress enter to continue, q enter to quit.> ");
        ans = chomp(readline())
        if ans == "q"
            break
        end
    end
end

function reviewdays(df::DataFrame)
    for it in eachrow(df)
        display(it)
        print("\nPress enter to continue, q enter to quit.> ");
        ans = chomp(readline())
        if ans == "q"
            break
        end
    end
end


#################################################################################
#  Epidemiological Stats -- very preliminary
#################################################################################

# outcomes per agegrp
function virus_outcome(series, locale; conds=[recovered, dead, totinfected, unexposed],
         agegrp=totalcol, base=:infected)  # denom in (:infected, :pop, :none)

    n = size(series[locale][:cum],1)
    outcomes = Dict{String, Float64}()  # TODO should we have integer outcomes for totals when base=:none?

    # each denominator for data summary
    total_pop = series[locale][:cum][1, map2series.unexposed[agegrp]] + series[locale][:cum][1, map2series.infectious[agegrp]]
    total_infected = series[locale][:cum][end, map2series.totinfected[agegrp]]

    denom = if base == :pop 
                total_pop 
            elseif base == :infected
                total_infected
            else  # :none or wrong entry
                1
            end

    for cond in conds
        outcomes[condnames[cond]] = series[locale][:cum][n, map2series[cond][agegrp]] / denom
    end

    return outcomes
end


function onecond(series, locale, cond; case=:new, agegrp=totalcol, filt=:pos)
    datacol = series[locale][case][:,map2series[cond][agegrp]]
    if filt == :pos
        datacol[datacol .> 0]
    else
        datacol
    end
end

###########################################################################################
#  Plotting
###########################################################################################


function cumplot(series, locale, plcols=[unexposed, infectious, recovered, dead]; 
    days="all", geo=[], thm=:wong2)

    pyplot()
    # theme(:ggplot2, foreground_color_border =:black, reuse = false)
    theme(thm, foreground_color_border=:black, 
          tickfontsize=9, gridlinewidth=1)

    !(typeof(plcols) <: Array) && (plcols = [plcols])

    # the data
    n = size(series[locale][:cum],1)
    days = days == "all" ? (1:n) : days
    cumseries = series[locale][:cum][days, [map2series[i][totalcol] for i in plcols]]

    # labels and annotations
    labels = [titlecase(condnames[i]) for i in plcols]
    labels = reshape([labels...], 1, length(labels))
    people = if !isempty(geo)
                geo[geo[:,fips] .== locale, popsize][1]
             else # this will off by a tiny bit because of rounding
                series[locale][:cum][1, map2series[unexposed][totalcol]] + series[locale][:cum][1,map2series[infectious][totalcol]]
             end   
    cityname = !isempty(geo) ? geo[geo[:,fips] .== locale, city][1] : ""
    died = series[locale][:cum][end, map2series[dead][totalcol]]
    # infected = series[locale][:cum][1,map2series[unexposed][totalcol]] - series[locale][:cum][end,map2series[unexposed][totalcol]]
    infected = people - series[locale][:cum][end,map2series[unexposed][totalcol]]
    recovered = infected - died
    unexp = people - infected

    firstseries = plcols[1]
    half_yscale = floor(Int, maximum(series[locale][:cum][:,map2series[firstseries][totalcol]]) * 0.7)
    co_pal = length(plcols) == 2 ? [theme_palette(thm)[2], theme_palette(thm)[4]] : theme_palette(thm)
 

    # the plot
    plot(   days, cumseries[days,1:end], 
            size = (700,500),
            label = labels, 
            lw=2.3,
            title = "Covid for $people people in $cityname over $n days\nActive Cases for Each Day",
            xlabel = "Simulation Days",
            ylabel = "People",
            legendfontsize = 10,
            color_palette = co_pal,
            reuse = false,
            annotate = ((6,half_yscale,
                Plots.text("Died: $died\nInfected: $infected\nRecovered: $recovered\nUnexposed: $unexp", 
                    11, :left)))
        )
    # annotate!((6,half_yscale,Plots.text("Died: $died\nInfected: $infected", 10, :left)))
    # gui()
end


function newplot(series, locale, plcols=[infectious]; days="all")

    # pyplot()
    theme(:ggplot2, foreground_color_border =:black)

    !(typeof(plcols) <: Array) && (plcols = [plcols])

    # the data and labels
    n = size(series[locale][:new],1)
    days = days == "all" ? (1:n) : days
    newseries = series[locale][:new][days, [map2series[i][totalcol] for i in plcols]]
    labels = [titlecase(condnames[i]) for i in plcols]
    labels = reshape([labels...], 1, length(labels))
    people = series[locale][:cum][1, map2series[unexposed][totalcol]] + series[locale][:cum][1,map2series[infectious][totalcol]]

    # the plot
    groupedbar( days, newseries[days,1:end], 
                size = (700,500),
                label = labels, 
                lw=0.2,
                bar_width=1,
                title = "Covid Daily Change for $people people over $n days",
                xlabel = "Simulation Days",
                yaxis = ("People"),
                reuse =false
        )
    gui()
end


function day2df(spreadq::Array)
    spreadseries = DataFrame(spreadq)

    spreadseries[!, :cuminfected] .= zeros(Int, size(spreadseries,1))
    spreadseries[1, :cuminfected] = copy(spreadseries[1,:infected])
    for i = 2:size(spreadseries,1)
       spreadseries[i,:cuminfected] = spreadseries[i-1,:cuminfected] + spreadseries[i,:infected]
    end

    return spreadseries
end


function dayplot(spreadq, plseries=[])
    dayplot(DataFrame(spreadq), plseries)
end


function dayplot(spreadseries::DataFrame, plseries=[])
    
    pyplot()
    theme(:ggplot2, foreground_color_border =:black)
    
    pl = bar(   spreadseries[!,:day], spreadseries[!,:infected],label="Infected", 
            lw=0.2,
            bar_width=1,
            size = (700,300),
            dpi=180,
            xlabel="Simulation Days", 
            ylabel="People", 
            title="Daily Spread of Covid",
            bg_legend=:white)
    
    for addlseries in plseries
        lbl = titlecase(string(addlseries))
        plot!(spreadseries[!,:day], spreadseries[!,addlseries],label=lbl, lw=2)
    end
    # gui()  # force instant plot window
    return pl
end


function day_animate2(spreadseries)
    n = size(spreadseries,1)
    # daymat = Matrix(spreadseries)

    xd = spreadseries[1:5,:]

    topy = max(maximum(spreadseries[!,:spreaders]),maximum(spreadseries[!,:contacts]),
                maximum(spreadseries[!,:touched]),maximum(spreadseries[!,:infected]) )

    @df xd plot(:day, [:spreaders :contacts :touched :infected], color=^([:red :blue :green :orange]),
                labels=^(["Spreaders" "Contacts" "Touched" "Infected"]),dpi=200, lw=2,ylim=(0,topy))

    for i = 5:2:n
        xd = spreadseries[i-2:i,:]

        @df xd plot!(:day, [:spreaders :contacts :touched :infected], color=^([:red :blue :green :orange]),
                 labels=false, dpi=200, lw=2, ylim=(0,3e4))
        gui()

        if i < round(Int, n/4)
            sleep(0.3)
        elseif i < round(Int,n/2)
            sleep(0.1)
        else
            sleep(.001)
        end
        # print("\nPress enter to continue, q enter to quit.> ");
        # ans = chomp(readline()) 
        # if ans == "q"
        #     break
        # end    
    end
end

function catplot()
    # groupedbar(datpct', bar_position=:stack,label=labels)
    # plot!(xticks=(1:3,["one", "two", "three"]))
    #= julia> datpct'
                3×5 LinearAlgebra.Adjoint{Float64,Array{Float64,2}}:
                 0.32036   0.30348    0.141939  0.155273  0.0789479
                 0.199046  0.413894   0.201274  0.054734  0.131053
                 0.252381  0.0677606  0.26775   0.161113  0.250996
    =#
end

# Plots.AnimatedGif("/var/folders/mf/73qj_8c91dzg4sw459_7mchm0000gn/T/jl_Js4px6.gif")