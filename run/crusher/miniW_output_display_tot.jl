############################################################
# Juilia miniWeather plotting script
############################################################

### USER DEFINE SECTION ####################################

file_name = "./output.nc"   # miniWeather output file
fig_dir_name = "./figures"  # Figure directory
make_ani = true             # Make animation

############################################################

using NCDatasets
using PyPlot
using Printf

XLEN = Float64(2.E4) # Miniweather X-dir length
ZLEN = Float64(1.E4) # Miniweather Z-dir height

extent=[0,XLEN/1000.0,0,ZLEN/1000.0]

# Open NetCDF output file with a create mode
ds = Dataset(file_name,"r")

ds_t = ds["t"]
ds_dens = ds["dens"]
ds_uwnd = ds["uwnd"]
ds_wwnd = ds["wwnd"]
ds_theta = ds["theta"]
ds_shum = ds["shumid"]

nt = ds.dim["t"]
nx = ds.dim["x"]
nz = ds.dim["z"]

# Make a figure directory (Remove & Make if exist)
if isdir(fig_dir_name)
   run(Cmd(`rm -rf $fig_dir_name`))
   run(Cmd(`mkdir -p $fig_dir_name`))
else
   run(Cmd(`mkdir -p $fig_dir_name`))
end

###############################
# Time loop
dens_max_t = zeros(Float64,nt)
uwnd_max_t = zeros(Float64,nt)
wwnd_max_t = zeros(Float64,nt)
rhot_max_t = zeros(Float64,nt)
shum_max_t = zeros(Float64,nt)

dens_min_t = zeros(Float64,nt)
uwnd_min_t = zeros(Float64,nt)
wwnd_min_t = zeros(Float64,nt)
rhot_min_t = zeros(Float64,nt)
shum_min_t = zeros(Float64,nt)

dens_max_t[1] = maximum((ds_dens[:,:,1]))
uwnd_max_t[1] = maximum((ds_uwnd[:,:,1]))
wwnd_max_t[1] = maximum((ds_wwnd[:,:,1]))
rhot_max_t[1] = maximum((ds_theta[:,:,1]))
shum_max_t[1] = maximum((ds_shum[:,:,1]))

dens_min_t[1] = minimum((ds_dens[:,:,1]))
uwnd_min_t[1] = minimum((ds_uwnd[:,:,1]))
wwnd_min_t[1] = minimum((ds_wwnd[:,:,1]))
rhot_min_t[1] = minimum((ds_theta[:,:,1]))
shum_min_t[1] = minimum((ds_shum[:,:,1]))

for n in 2:nt
    dens_max_t[n] = max(dens_max_t[n-1],maximum((ds_dens[:,:,n])))
    uwnd_max_t[n] = max(uwnd_max_t[n-1],maximum((ds_uwnd[:,:,n])))
    wwnd_max_t[n] = max(wwnd_max_t[n-1],maximum((ds_wwnd[:,:,n])))
    rhot_max_t[n] = max(rhot_max_t[n-1],maximum((ds_theta[:,:,n])))
    shum_max_t[n] = max(shum_max_t[n-1],maximum((ds_shum[:,:,n])))

    dens_min_t[n] = min(dens_min_t[n-1],minimum((ds_dens[:,:,n])))
    uwnd_min_t[n] = min(uwnd_min_t[n-1],minimum((ds_uwnd[:,:,n])))
    wwnd_min_t[n] = min(wwnd_min_t[n-1],minimum((ds_wwnd[:,:,n])))
    rhot_min_t[n] = min(rhot_min_t[n-1],minimum((ds_theta[:,:,n])))
    shum_min_t[n] = min(shum_min_t[n-1],minimum((ds_shum[:,:,n])))
end

dens_max = (maximum(dens_max_t)) * 0.5
uwnd_max = (maximum(uwnd_max_t)) * 0.5
wwnd_max = (maximum(wwnd_max_t)) * 0.5
rhot_max = (maximum(rhot_max_t)) * 0.5
shum_max = (maximum(shum_max_t)) * 0.5

dens_min = (minimum(dens_min_t)) * 0.5
uwnd_min = (minimum(uwnd_min_t)) * 0.5
wwnd_min = (minimum(wwnd_min_t)) * 0.5
rhot_min = (minimum(rhot_min_t)) * 0.5
shum_min = (minimum(shum_min_t)) * 0.5

###############################

# Time loop
for n in 1:nt
    println("Processing time: ",n," / ",nt)
 
    time = ds_t[n]
    dens = ds_dens[:,:,n]
    uwnd = ds_uwnd[:,:,n]
    wwnd = ds_wwnd[:,:,n]
    rhot = ds_theta[:,:,n]
    shum = ds_shum[:,:,n]

    it = @sprintf("%.4i",n)
    ctime = @sprintf("%6.1f",time)

    # Display
    fig = figure("out",figsize=(13,5.5))
    subplots_adjust(hspace=0.1)
    suptitle(string("Julia Miniweather
            \nnx=",nx,"; nz=",nz,"; Time=",ctime," s"))
    subplot(231)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Density")
       imshow(transpose(reverse(dens,dims=2)),extent=extent,cmap="jet",vmin=dens_min,vmax=dens_max)
    subplot(232)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Potential temperature")
       imshow(transpose(reverse(rhot,dims=2)),extent=extent,cmap="jet",vmin=rhot_min,vmax=rhot_max)
    subplot(233)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Specific humidity")
       imshow(transpose(reverse(shum,dims=2)),extent=extent,cmap="jet",vmin=shum_min,vmax=shum_max)
    subplot(234)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("U-wind component")
       imshow(transpose(reverse(uwnd,dims=2)),extent=extent,cmap="bwr",vmin=uwnd_min,vmax=uwnd_max)
    subplot(235)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("W-wind component")
       imshow(transpose(reverse(wwnd,dims=2)),extent=extent,cmap="bwr",vmin=wwnd_min,vmax=wwnd_max)

    fig_name = string(fig_dir_name,"/fig_",it,".png")
    savefig(fig_name,bbox_inches="tight",dpi=120)
    close(fig)
end

# Make figures as animation
if make_ani
   run(Cmd(`convert -delay 0 -loop 0 $fig_dir_name/\*.png anim.gif`))
end
