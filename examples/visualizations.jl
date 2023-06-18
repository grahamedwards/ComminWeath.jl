using ComminWeath, CairoMakie

 
const cwcolorscheme = [:hotpink2,:midnightblue, :lightslategray,:coral, :cyan3,:plum3]

function plotUseries(a::Vector,t::Vector;f=Figure(), cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind(),meas::NamedTuple=(;))
    (@isdefined set_theme!) || error("You must load Makie for this function to work, e.g. using CairoMakie")
    set_theme!(; palette=(; color=cwcolorscheme))
    isempty(cU) || @assert length(cU)==length(a) "cU and a must be the same length"

    A234 = zeros(Float64, length(a),length(t))
    A230 = copy(A234)
    cUout = copy(A234)
    cwouts = Vector{NamedTuple}(undef,length(a))

    ax = Axis(f[1,1], xlabel="(²³⁰Th/²³⁸U)",ylabel="(²³⁴U/²³⁸U)",xgridvisible=false,ygridvisible=false)
    if !isempty(meas)
        errorbars!(ax, meas.r08, meas.r48, meas.u08, whiskerwidth = 0, direction = :x,linewidth=2,color=:black)
        errorbars!(ax, meas.r08, meas.r48, meas.u48, whiskerwidth = 0, direction = :y,linewidth=2,color=:black)
    end
    for i in eachindex(a)
        isempty(cU) || (d.cU=cU[i])
        cwouts[i] = comminweath(a[i],g,d,wx,r)
        for j in eachindex(t)
            Uout = drawdate(t[j],cwouts[i])
            A234[i,j] = Uout.A234
            A230[i,j] = Uout.A230
            cUout[i,j] = Uout.cU
        end
        tjmax = searchsortedfirst(cwouts[i].t,maximum(t)*1e3)
        lines!(ax,cwouts[i].A230[1:tjmax],cwouts[i].A234[1:tjmax],linewidth=1,color=:gray60,)
    end

    for i in eachindex(t)
        x = view(A230,:,i)
        y = view(A234,:,i)
        if !iszero(t[i])
            reg = linreg(x,y)
            lines!(ax,[first(x),last(x)],[reg.b+reg.m*first(x),reg.b+reg.m*last(x)],label=string(Int(t[i])," ka"),linewidth=2)
        end
        scatter!(ax,x,y,markersize=10,color=:gray40)
    end
    Legend(f[1,1],ax,tellheight=false,tellwidth=false,halign=:left,valign=:top,margin=(10,10,10,10))
    f
end

function plotslopes(t::AbstractVector;f=Figure(), a::Vector=[10.,20.,30.,40.],cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind())
    (@isdefined set_theme!) || error("You must load Makie for this function to work, e.g. using CairoMakie")
    slopes = calcslopes(t,a=a,cU=cU,g=g,d=d,wx=wx,r=r)
    Axis(f[1,1], xlabel="Age (ka)", ylabel="slope",xgridvisible=false,ygridvisible=false)
    lines!(t,slopes,color=:black,linewidth=2)
    f
end

function plotcU(a::Vector,t::Vector;f=Figure(), cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind(), meas::Vector=[])
    (@isdefined set_theme!) || error("You must load Makie for this function to work, e.g. using CairoMakie")
    set_theme!(; palette=(; color=[:gray40, cwcolorscheme...]))
    ax=Axis(f[1,1], xlabel="Grain diameter (μm)",ylabel="[U] (μg g⁻¹)",xgridvisible=false,ygridvisible=false)
    cUout = calcU(a,t,cU=cU, g=g, d=d,wx=wx,r=r)
    cUout .*= 1e-3

    if !isempty(meas)
        @assert length(a) == length(meas)
        for i in eachindex(a)
            d = a[i]
            upper=meas[i][2]
            lower=meas[i][1]
            band!(ax,[d-1,d+1], [lower,lower], [upper,upper],color=(:lightsteelblue,0.4))
        end
    end

    for i in eachindex(t)
        scatterlines!(ax,a,view(cUout,:,i),label=string(Int(t[i])," ka"))
    end
    Legend(f[1,1],ax,tellheight=false,tellwidth=false,halign=:left,valign=:top,margin=(10,10,10,10))
    f
end

function plotauthreplace(a::Vector,t::Vector;f=Figure(), cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind(), meas::Vector=[])
    (@isdefined set_theme!) || error("You must load Makie for this function to work, e.g. using CairoMakie")
    set_theme!(; palette=(; color=[:gray40, cwcolorscheme...]))
    ax=Axis(f[1,1], xlabel="Grain diameter (μm)",ylabel="Authigenic replacement (%)",xgridvisible=false,ygridvisible=false)
    cUout = calcU(a,t,cU=cU, g=g, d=d,wx=wx,r=r)
    cUₒ= ifelse(isempty(cU),fill(d.cU,length(a)),cU)

    if !isempty(meas)
        @assert length(a) == length(meas)
        for i in eachindex(a)
            d = a[i]
            upper=meas[i][2]
            lower=meas[i][1]
            band!(ax,[d-1,d+1], [lower,lower], [upper,upper],color=(:lightsteelblue,0.4))
        end
    end

    for i in eachindex(t)
        authpct =  100 .* ( view(cUout,:,i) .- cUₒ ) ./ (wx.cU .- cUₒ)
        scatterlines!(ax,a,authpct,label=string(Int(t[i])," ka"))
    end
    Legend(f[1,1],ax,tellheight=false,tellwidth=false,halign=:right,valign=:top,margin=(10,10,10,10))
    f
end