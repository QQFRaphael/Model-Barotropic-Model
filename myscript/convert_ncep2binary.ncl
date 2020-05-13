; This script is used to rewrite NCEP-DOE reanalysis II
; to binary format. The output of this script will be used
; to create background flow of the barotropical model
; Please note only 30 years can be used, because
; code/GetVorticity.f only read 30 years if not changed
; if more years is needed, modifiy Line12 accordingly

load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/calendar_decode2.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

yearstart = 1979
yearend = 2018
nyear = yearend - yearstart + 1

fout1 = "umon200.dat"
fout2 = "vmon200.dat"

f1 = addfile("../data/uwnd.mon.mean.nc","r")
f2 = addfile("../data/vwnd.mon.mean.nc","r")

year=calendar_decode2(f1->time,0)
year_idx=ind(year(:,0).ge.(yearstart).and.year(:,0).le.(yearend))
level = f1->level
lev_idx = ind(level.eq.200.0)

u=short2flt(f1->uwnd(year_idx,lev_idx,::-1,:))
v=short2flt(f2->vwnd(year_idx,lev_idx,::-1,:))

system("rm -f " + fout1)
system("rm -f " + fout2)

fbindirwrite(fout1,u)
fbindirwrite(fout2,v)
