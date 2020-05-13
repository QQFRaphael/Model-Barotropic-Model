; This script is used to plot model results.
; As this is a barotropical model, only winds and HGT can be shown

load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

case   = "test"
starts = 30
ends   = 39

f_clim = addfile("../output/"+case+"/stream_ctl.nc", "r")
f_frc  = addfile("../output/"+case+"/stream_frc.nc", "r")

sf_clim = f_clim->sf
sf_frc  = f_frc->sf
lat = f_clim->lat

fg = doubletofloat(2.0*5.29*10^-5*sin(lat*3.14159/180.)/9.8)

diff = sf_frc - sf_clim
diff = diff - conform(diff, dim_avg_n_Wrap(diff, 2), (/0,1/))
diff = diff * conform(diff, fg, 1)
copy_VarMeta(sf_clim, diff)

wks = gsn_open_wks("png",case)

res                              = True

res@gsnDraw                      = False
res@gsnFrame                     = False
res@gsnAddCyclic                 = True

res@gsnLeftString                = ""  
res@gsnRightString               = ""
res@tiXAxisString                = ""
res@tiYAxisString                = ""

res@cnFillOn                     = False
res@cnLinesOn                    = True 
res@cnInfoLabelOn                = False
res@cnLineLabelsOn               = True 
res@cnLevelSelectionMode         = "ExplicitLevels" 
res@cnLevels                     = ispan(-700,700,50)
res@cnLineLabelPlacementMode     = "constant"    
res@cnLineLabelInterval          =  2 
res@gsnContourNegLineDashPattern = 1
res@cnLineThicknessF             = 2
res@cnLineLabelFormat            = "@*^sg" 
res@cnSmoothingOn                = True 
res@cnLineLabelFontHeightF       = 0.011 

res@mpCenterLonF                 = 60
res@mpMinLatF                    = 0  
res@mpMaxLatF                    = 90   
res@mpMinLonF                    = -110   
res@mpMaxLonF                    = 220 

plot = gsn_csm_contour_map(wks,dim_avg_n_Wrap(diff(starts:ends,:,:), 0),res) 

draw(wks)
frame(wks)
