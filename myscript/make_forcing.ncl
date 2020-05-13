; This script is used to make forcing fields

load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

cent_lat = 37.5
cent_lon = 345.0
len_lon = 15.0
len_lat = 15.0

div_intensity = 7e-6
vor_intensity = 0.0

month_idx = (/5,6,7/)

f = addfile("../data/Gaussian_R40.nc", "r")

div = f->div
vor = f->vor
vpt = f->vpt
str = f->str

lat = f->LAT
lon = f->LON

; the shape, intensity of forcing should be created by yourself

div_mean = dim_avg_n_Wrap(div(month_idx,:,:), 0)
vor_mean = dim_avg_n_Wrap(vor(month_idx,:,:), 0)

shape = new((/dimsizes(div_mean)/), "float")
lat2d = conform(div_mean, lat, 0)
lon2d = conform(div_mean, lon, 1)

shape = exp(-(lat2d-cent_lat)*(lat2d-cent_lat)/len_lat)*exp(-(lon2d-cent_lon)*(lon2d-cent_lon)/len_lon)  ;smth9(exp(-(lat2d-cent_lat)*(lat2d-cent_lat)/len_lat)*exp(-(lon2d-cent_lon)*(lon2d-cent_lon)/len_lon),0.5,-0.25,True)
shape = where((shape-0.01).ge.0.0, shape, 0.0)*10.0
copy_VarMeta(div_mean, shape)

div_frc = div_mean + shape*(div_intensity)
copy_VarMeta(div_mean, div_frc)

vor_frc = vor_mean + shape*(vor_intensity)
copy_VarMeta(vor_mean, vor_frc)

system("rm -rf ../data/frc_bs.nc")
out = addfile("../data/frc_bs.nc", "c")
out->div_mean = div_mean
out->div_frc = div_frc
out->vor_mean = vor_mean
out->vor_frc = vor_frc

fbindirwrite("../data/div.grads.clim",div_mean)
fbindirwrite("../data/div.grads.frc",div_frc)
fbindirwrite("../data/vor.grads.clim",vor_mean)
fbindirwrite("../data/vor.grads.frc",vor_frc)
