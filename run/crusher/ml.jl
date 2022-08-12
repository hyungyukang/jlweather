############################################################
# Juilia miniWeather plotting script
############################################################

### USER DEFINE SECTION ####################################

file_name = "./output.nc"   # miniWeather output file
fig_dir_name = "./figures"  # Figure directory
make_ani = false            # Make animation

############################################################

using NCDatasets
using PyPlot
using Printf
using Flux

XLEN = Float64(2.E4) # Miniweather X-dir length
ZLEN = Float64(1.E4) # Miniweather Z-dir height

extent=[0,XLEN/1000.0,0,ZLEN/1000.0]

# Open NetCDF output file with a create mode
ds = Dataset(file_name,"r")

ds_t = ds["t"]
ds_dens_bef =    ds["dens_bef"]
ds_uwnd_bef =    ds["uwnd_bef"]
ds_wwnd_bef =    ds["wwnd_bef"]
ds_theta_bef =  ds["theta_bef"]
ds_shum_bef =  ds["shumid_bef"]
ds_dens_aft =    ds["dens_aft"]
ds_uwnd_aft =    ds["uwnd_aft"]
ds_wwnd_aft =    ds["wwnd_aft"]
ds_theta_aft =  ds["theta_aft"]
ds_shum_aft =  ds["shumid_aft"]

nt = ds.dim["t"]
nx = ds.dim["x"]
nz = ds.dim["z"]

# Make a figure directory (Remove & Make if exist)
#if isdir(fig_dir_name)
#   run(Cmd(`rm -rf $fig_dir_name`))
#   run(Cmd(`mkdir -p $fig_dir_name`))
#else
#   run(Cmd(`mkdir -p $fig_dir_name`))
#end


x_vec = zero(Float64,nz*5)
y_vec = zero(Float64,nz*5)

# Time loop
for n in 1:nt
    println("Processing time: ",n," / ",nt)
 
    time = ds_t[n]
    dens_bef =  ds_dens_bef[:,:,n]
    uwnd_bef =  ds_uwnd_bef[:,:,n]
    wwnd_bef =  ds_wwnd_bef[:,:,n]
    rhot_bef = ds_theta_bef[:,:,n]
    shum_bef =  ds_shum_bef[:,:,n]

    dens_aft =  ds_dens_aft[:,:,n]
    uwnd_aft =  ds_uwnd_aft[:,:,n]
    wwnd_aft =  ds_wwnd_aft[:,:,n]
    rhot_aft = ds_theta_aft[:,:,n]
    shum_aft =  ds_shum_aft[:,:,n]



#   it = @sprintf("%.4i",n)
#   ctime = @sprintf("%6.1f",time)

#   # Display
#   fig = figure("out",figsize=(13,5.5))
#   subplots_adjust(hspace=0.1)
#   suptitle(string("Julia Miniweather
#           \nnx=",nx,"; nz=",nz,"; Time=",ctime," s"))
#   subplot(231)
#      xlabel("X (km)")
#      ylabel("Z (km)")
#      title("Density")
#      imshow(transpose(reverse(dens,dims=2)),extent=extent,cmap="jet")
#   subplot(232)
#      xlabel("X (km)")
#      ylabel("Z (km)")
#      title("Potential temperature")
#      imshow(transpose(reverse(rhot,dims=2)),extent=extent,cmap="jet")
#   subplot(233)
#      xlabel("X (km)")
#      ylabel("Z (km)")
#      title("Specific humidity")
#      imshow(transpose(reverse(shum,dims=2)),extent=extent,cmap="jet")
#   subplot(234)
#      xlabel("X (km)")
#      ylabel("Z (km)")
#      title("U-wind component")
#      imshow(transpose(reverse(uwnd,dims=2)),extent=extent,cmap="bwr")
#   subplot(235)
#      xlabel("X (km)")
#      ylabel("Z (km)")
#      title("W-wind component")
#      imshow(transpose(reverse(wwnd,dims=2)),extent=extent,cmap="bwr")

#   fig_name = string(fig_dir_name,"/fig_",it,".png")
#   savefig(fig_name,bbox_inches="tight",dpi=120)
#   close(fig)
end

# Make figures as animation
if make_ani
   run(Cmd(`convert -delay 0 -loop 0 $fig_dir_name/\*.png anim.gif`))
end
