############################################################
# Juilia miniWeather plotting script
############################################################

### USER DEFINE SECTION ####################################

file_name_dry = "./output_dry.nc"   # miniWeather output file
file_name_moist_np = "./output_moist_noPrecp.nc"   # miniWeather output file
file_name_moist_p = "./output_moist_Precp.nc"   # miniWeather output file
fig_dir_name = "./figures"  # Figure directory
make_ani = true            # Make animation

############################################################

using NCDatasets
using PyPlot
using Printf

XLEN = Float64(2.E4) # Miniweather X-dir length
ZLEN = Float64(1.E4) # Miniweather Z-dir height

extent=[0,XLEN/1000.0,0,ZLEN/1000.0]

# Open NetCDF output file with a create mode
ds_d = Dataset(file_name_dry,"r")
ds_m_np = Dataset(file_name_moist_np,"r")
ds_m_p = Dataset(file_name_moist_p,"r")

ds_t = ds_d["t"]
ds_dens_d = ds_d["dens"]
ds_uwnd_d = ds_d["uwnd"]
ds_wwnd_d = ds_d["wwnd"]
ds_thet_d = ds_d["theta"]

ds_dens_m_np = ds_m_np["dens_aft"]
ds_uwnd_m_np = ds_m_np["uwnd_aft"]
ds_wwnd_m_np = ds_m_np["wwnd_aft"]
ds_thet_m_np = ds_m_np["theta_aft"]
ds_shum_m_np = ds_m_np["shumid_aft"]

ds_dens_m_p = ds_m_p["dens_aft"]
ds_uwnd_m_p = ds_m_p["uwnd_aft"]
ds_wwnd_m_p = ds_m_p["wwnd_aft"]
ds_thet_m_p = ds_m_p["theta_aft"]
ds_shum_m_p = ds_m_p["shumid_aft"]


nt = ds_d.dim["t"]
nx = ds_d.dim["x"]
nz = ds_d.dim["z"]

# Make a figure directory (Remove & Make if exist)
if isdir(fig_dir_name)
   run(Cmd(`rm -rf $fig_dir_name`))
   run(Cmd(`mkdir -p $fig_dir_name`))
else
   run(Cmd(`mkdir -p $fig_dir_name`))
end

pvmax=5.0
pvmin=0.0
qvmax=0.005
qvmin=0.0

# Time loop
for n in nt-1:nt
    println("Processing time: ",n," / ",nt)
 
    time = ds_t[n]

    dens_d = ds_dens_d[:,:,n]
    uwnd_d = ds_uwnd_d[:,:,n]
    wwnd_d = ds_wwnd_d[:,:,n]
    rhot_d = ds_thet_d[:,:,n]

    dens_m_np = ds_dens_m_np[:,:,n]
    uwnd_m_np = ds_uwnd_m_np[:,:,n]
    wwnd_m_np = ds_wwnd_m_np[:,:,n]
    rhot_m_np = ds_thet_m_np[:,:,n]
    shum_m_np = ds_shum_m_np[:,:,n]

    dens_m_p = ds_dens_m_p[:,:,n]
    uwnd_m_p = ds_uwnd_m_p[:,:,n]
    wwnd_m_p = ds_wwnd_m_p[:,:,n]
    rhot_m_p = ds_thet_m_p[:,:,n]
    shum_m_p = ds_shum_m_p[:,:,n]

    it = @sprintf("%.4i",n)
    ctime = @sprintf("%6.1f",time)

    # Display
    fig = figure("out",figsize=(13,8.0))
    subplots_adjust(hspace=0.1)
    suptitle(string("Julia Miniweather
            \nnx=",nx,"; nz=",nz,"; Time=",ctime," s"))
    subplot(331)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Potential temperature")
       imshow(transpose(reverse(rhot_d,dims=2)),extent=extent,cmap="jet",vmin=pvmin,vmax=pvmax)
    subplot(332)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Potential temperature")
       imshow(transpose(reverse(rhot_m_np,dims=2)),extent=extent,cmap="jet",vmin=pvmin,vmax=pvmax)
    subplot(333)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Potential temperature")
       imshow(transpose(reverse(rhot_m_p,dims=2)),extent=extent,cmap="jet",vmin=pvmin,vmax=pvmax)
    subplot(334)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Specific humidity")
       imshow(transpose(reverse(rhot_d*0,dims=2)),extent=extent,cmap="jet",vmin=qvmin,vmax=qvmax)
       text(0.5,0.5,"Const.=0",color="white")
    subplot(335)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Specific humidity")
       imshow(transpose(reverse(shum_m_np,dims=2)),extent=extent,cmap="jet",vmin=qvmin,vmax=qvmax)
    subplot(336)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Specific humidity")
       imshow(transpose(reverse(shum_m_p,dims=2)),extent=extent,cmap="jet",vmin=qvmin,vmax=qvmax)
    subplot(337)
       xlabel("X (km)")
       ylabel("Z (km)")
       title("Specific humidity")
       quiver((transpose(reverse(uwnd_m_p,dims=2))),(transpose(reverse(wwnd_m_p,dims=2))))

    fig_name = string(fig_dir_name,"/fig_",it,".png")
    savefig(fig_name,bbox_inches="tight",dpi=120)
    close(fig)
end

# Make figures as animation
if make_ani
   run(Cmd(`convert -delay 0 -loop 0 $fig_dir_name/\*.png anim.gif`))
end
