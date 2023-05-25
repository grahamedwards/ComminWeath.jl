using CairoMakie, ColorSchemes

function plotUseries(a::Vector,t::Vector;f=Figure(), cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind(),meas::NamedTuple=(;))
    set_theme!(; palette=(; color=[:hotpink2,:midnightblue,ColorSchemes.seaborn_colorblind6...]))
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


function calcslopes(t::AbstractVector; a::Vector=[10.,20.,30.,40.],cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind())
    A234 = zeros(Float64, length(a),length(t))
    A230 = copy(A234)
    cwouts = Vector{NamedTuple}(undef,length(a))
    slopes = similar(t,Float64)

    for i in eachindex(a)
        isempty(cU) || (d.cU=cU[i])
        cwouts[i] = comminweath(a[i],g,d,wx,r)
        for j in eachindex(t)
            Uout = drawdate(t[j],cwouts[i])
            A234[i,j] = Uout.A234
            A230[i,j] = Uout.A230
        end
    end
    for i in eachindex(t)
        x = view(A230,:,i)
        y = view(A234,:,i)
        if !iszero(t[i])
            slopes[i] = linreg(x,y).m
        end
    end
    slopes
end

function plotslopes(t::AbstractVector;f=Figure(), a::Vector=[10.,20.,30.,40.],cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind())
    slopes = calcslopes(t,a=a,cU=cU,g=g,d=d,wx=wx,r=r)
    Axis(f[1,1], xlabel="Age (ka)", ylabel="slope",xgridvisible=false,ygridvisible=false)
    lines!(t,slopes,color=:black,linewidth=2)
    f
end


function calcU(a::Vector,t::Vector;cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind())
    cUout = zeros(Float64, length(a),length(t))
    cwouts = Vector{NamedTuple}(undef,length(a))
    for i in eachindex(a)
        isempty(cU) || (d.cU=cU[i])
        cwouts[i] = comminweath(a[i],g,d,wx,r)
        for j in eachindex(t)
            cUout[i,j] = drawdate(t[j],cwouts[i]).cU
        end
    end
    cUout
end

function plotcU(a::Vector,t::Vector;f=Figure(), cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind())
    set_theme!(; palette=(; color=[:gray40,:hotpink2,:midnightblue,ColorSchemes.seaborn_colorblind6...]))
    ax=Axis(f[1,1], xlabel="Grain diameter (μm)",ylabel="[U] (μg g⁻¹)",xgridvisible=false,ygridvisible=false)
    cUout = calcU(a,t,cU=cU, g=g, d=d,wx=wx,r=r)
    cUout .*= 1e-3
    for i in eachindex(t)
        scatterlines!(ax,a,view(cUout,:,i),label=string(Int(t[i])," ka"))
    end
    Legend(f[1,1],ax,tellheight=false,tellwidth=false,halign=:left,valign=:top,margin=(10,10,10,10))
    f
end


function plotauthreplace(a::Vector,t::Vector;f=Figure(), cU::Vector=[],g::Grain=Grain(),d::Detrital=Detrital(),wx::WxAuth=WxAuth(),r::Rind=Rind())
    set_theme!(; palette=(; color=[:gray40,:hotpink2,:midnightblue,ColorSchemes.seaborn_colorblind6...]))
    ax=Axis(f[1,1], xlabel="Grain diameter (μm)",ylabel="[U] (μg g⁻¹)",xgridvisible=false,ygridvisible=false)
    cUout = calcU(a,t,cU=cU, g=g, d=d,wx=wx,r=r)
    println(cUout)
    cUₒ= ifelse(isempty(cU),fill(d.cU,length(a)),cU)
    for i in eachindex(t)
        authpct =  100 .* ( view(cUout,:,i) .- cUₒ ) ./ (wx.cU .- cUₒ)
        scatterlines!(ax,a,authpct,label=string(Int(t[i])," ka"))
    end
    Legend(f[1,1],ax,tellheight=false,tellwidth=false,halign=:right,valign=:top,margin=(10,10,10,10))
    f
end