  # Check stabilities
  stable=true
  radius=Hill_radius(avg[12],(avg[11].*CGS.MSUN/CGS.MEARTH),calc_ecc(avg[14],avg[15]),CGS.MSUN)
  ab_mutual_radius=mutual_Hill(avg[2],avg[1].*CGS.MSUN/CGS.MEARTH,CGS.MSUN,avg[7],avg[6].*CGS.MSUN/CGS.MEARTH)
  bc_mutual_radius=mutual_Hill(avg[7],avg[6].*CGS.MSUN/CGS.MEARTH,CGS.MSUN,avg[12],avg[11].*CGS.MSUN/CGS.MEARTH)
  println("Planet a-b mutual Hill_radius: ",ab_mutual_radius)
  println("Planet b-c mutual Hill_radius: ",bc_mutual_radius)
  println("Planet c Hill radius: ",radius)
  # period_ratios=[]
  # periods=[avg[2],avg[7],avg[12]]
  # for i=1:nplanet
  #   if i<3
  #     append!(period_ratios,periods[i+1]/periods[i])
  #   end
  # end
nbins=10
N=10000
x,y=randn((N),2)
ax=figure()
xyhist,xbins,ybins=ax.hist2d(x,y,bins=nbins)
xbin_edges,ybin_edges=np.meshgrid(collect(range(minimum(x),maximum(x),length=nbins)),collect(range(minimum(y),maximum(y),length=nbins)))
ax.contour(xbin_edges,ybin_edges,xyhist)

## make corner plot to compare posterior distributions of 2 planets 
  subplots_adjust(hspace=0.05,wspace=0.05)
  N=length(xs)
  for col_item in 1:N,row_item in 1:N
    if col_item==row_item
      fig.add_subplot(N,N,col_item)
      axvline(true_xs[col_item],linestyle="--",color="black")
      hist(xs[col_item],bins=nbins,histtype="step",color="black")
      axvspan(quantile(xs[col_item],0.1587),quantile(xs[col_item],0.8413),color="limegreen",alpha=0.35)
    end

  grid_type_nplanet="widep4";sigma=30;nyear=30;
  mcfile=string("MCMC/",grid_type_nplanet,"_mcmc",sigma,"s",nyear,"yrs.jld2")
  m=jldopen(String(mcfile),"r")
  par_mcmc= m["par_mcmc"]
  lprob_mcmc=m["lprob_mcmc"]
  nwalkers=m["nwalkers"]
  nsteps=m["nsteps"]
  iburn=m["iburn"]
  indepsamples=m["indepsamples"]
  truem1=0.815
  trueec1=calc_evec1(0.00677323,131.53298)
  truees1=calc_evec2(0.00677323,131.53298)
  truep1=224.7007992#.-offset
  truee1=0.00677323
  truem2=1
  trueec2=calc_evec1(0.01671022,102.94719)
  truees2=calc_evec2(0.01671022,102.94719)
  truep2=365.2564#-offset #365.256355
  truee2=0.01671022
  truem3=0.1074
  trueec3=calc_evec1(0.09341233,336.04084)
  truees3=calc_evec2(0.09341233,336.04084)
  truep3=686.9795859
  truee3=0.09341233
  corner(par_mcmc[:,iburn:end,:],pname,50,[truem1,truep1,trueec1,truees1,truem2,truep2,trueec2,truees2])
function corner(params,labels,nbins,true_vals)
N=length(params)
ncols=N
nrows=N
columns=labels
fig=figure(figsize=(10,10))
# ax=plt.axis()
for row_ind=1:nrows, col_ind=1:ncols
    # ax.add_subplot()
  if col_ind==row_ind
    #make histograms
    for iparam=1:nparam
      # new_par_ind = (nparam%3)
    #   if iparam==new_par_ind
    #     break  
    #   else
    #     println(pname[iparam])
    #   end
    # end
      subplot(col_ind,row_ind,iparam)
      ax=gca()
      ax.hist((params[:,:,iparam]),bins=nbins,histtype="step",orientation="horizontal")
      # ax.axvspan(quantile(params[iparam],0.1587),quantile(params[iparam],0.8413),color="limegreen",alpha=0.35)
      # ax.axvline(true_vals[iparam],linestyle="--",color="black")
    end
  end

end  
end