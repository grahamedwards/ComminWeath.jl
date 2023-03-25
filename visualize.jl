## Plot Isochrons & Age Calculations
using Plots; gr() #load Plots with gr backend.
p1 = plot(A230,A234,
                     legend=false, label = Int.(tcom'/1000),
                     legendtitle="Age (ka)", marker = (:hexagon, 2, 0.6, :black))
xaxis!("(230Th/238U)")
ylabel!("(234U/238U)")
# Plot Evolutions
    for i=1:length(a)
        global p1 =plot!(A230_ev[i,length(tcom),:],A234_ev[i,length(tcom),:],
        color= :gray, linewidth = 0.1, label="")
    end
p1=plot!(meas.r08,meas.r48,xerr=meas.u08, yerr=meas.u48)
# Plot slope vs. age
p2 = plot(tcom/1000,slope_int[:,2], label ="") #, ylim=(0,0.6))
ylabel!("slope");xlabel!("Age (ka)")

# Plot [U]_detrital vs. grain size
p3 = plot(repeat(a,1,length(tcom))*1e6,cUd/1e3,
legend=false, label = Int.(tcom'/1000),
legendtitle="Age (ka)", marker = (:hexagon, 2, 0.6, :black))
ylabel!("[U] (ppm)");xlabel!("Grain size (um)")

# Calculate fraction clay from cUd
auth_pct = 100*(cUd .- repeat(U_dtr,1,length(tcom))) ./ (U_auth .- U_dtr)
p4 = plot(repeat(a,1,length(tcom))*1e6,auth_pct,
legend=false, label = Int.(tcom'/1000),
legendtitle="Age (ka)",marker = (:hexagon, 2, 0.6, :black))
ylabel!("Percent authigenic phase (%)");xlabel!("Grain size (um)")

plot_out=plot(p1,p2,p3,p4, layout=(2,2))
current()
savefig(string(fileout_prefix,run_name,".pdf"))
display(plot_out)