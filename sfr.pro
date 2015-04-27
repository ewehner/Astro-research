; Code written by Elizabeth Wehner in IDL
; Code loads magnitudes for galaxy NGC 3310 and calculates colors and luminosities
; Magnitudes were measured by E. Wehner using IRAF.
; 
; Final publication resulting from the following code can be found here:
; http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?2006MNRAS.371.1047W&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf
;



pro sfr

set_plot, 'ps'
device, filename='nice_plots.ps'

$ read in my magnitudes and calculate my colors

readcol, '~/data/n3310/final/all_mags_feb14.dat', u_ew, b_ew, v_ew, r_ew

UmB_ew = u_ew - b_ew
BmV_ew = b_ew - v_ew
BmR_ew = b_ew - r_ew
VmR_ew = v_ew - r_ew

$ read magnitudes and ages for z = z2 SD, used to represent the underlying galaxy (ug)

readcol, 'galev/sd_z2_salpeter_EL_fvm05/johnson_mag.dat', num, age, uj_ug, bj_ug, vj_ug, rj_ug, ij_ug, jj_ug, hj_ug, kj_ug, rc_ug, ic_ug
readcol, 'galev//sd_z2_salpeter_EL_fvm05/stat.dat', age, mass_gas, mass_stars, sfr

mass_tot_ug = mass_gas(2999) + mass_stars(2999)
print, "mass_tot_ug = ", mass_tot_ug

$ array element where t ~ 12Gyr, colors at this age are used to represent the underlying galaxy
step_ug = 333

$ Calculate the total luminosity based on magnitude, and then divide by the total mass in stars

Luj_ug = 3.86e26*10^(-0.4*(uj_ug(step_ug) - 4.76)) / mass_stars(2999)
Lbj_ug = 3.86e26*10^(-0.4*(bj_ug(step_ug) - 4.76)) / mass_stars(2999)
Lvj_ug = 3.86e26*10^(-0.4*(vj_ug(step_ug) - 4.76)) / mass_stars(2999)
Lrc_ug = 3.86e26*10^(-0.4*(rc_ug(step_ug) - 4.76)) / mass_stars(2999)

$ read magnitudes and ages for z = z0 E, used to represent the underlying galaxy (ug)
$
$readcol, 'galev/e_z0_salpeter_EL_fvm05/johnson_mag.dat', num, age_e, uj_e, bj_e, vj_e, rj_e, ij_e, jj_e, hj_e, kj_e, rc_e, ic_e
$UmB_e = uj_e - bj_e
$BmV_e = bj_e - vj_e
$VmR_e = vj_e - rc_e

$ read magnitudes and ages for the SSP used to represent the overlying burst

readcol, 'galev/ssp_z2_salpeter_EL_fvm05/johnson_mag.dat', num, age, uj_ssp, bj_ssp, vj_ssp, rj_ssp, ij_ssp, jj_ssp, hj_ssp, kj_ssp, rc_ssp, ic_ssp
readcol, 'galev/ssp_z2_salpeter_EL_fvm05/stat.dat', age, mass_gas_ssp, mass_stars_ssp, sfr_ssp

mass_tot_ssp = mass_gas_ssp(2999) + mass_stars_ssp(2999)
print, "mass_tot_ssp = ", mass_tot_ssp

sz = size(age)
num_steps = sz(1)

$ Calculate the total luminosities, based on the absolute magnitudes, and then divide by the stellar mass to normalize L to per solar mass

Luj_ssp = 3.86e26*10^(-0.4*(uj_ssp - 4.76)) / mass_stars_ssp(2999)
Lbj_ssp = 3.86e26*10^(-0.4*(bj_ssp - 4.76)) / mass_stars_ssp(2999)
Lvj_ssp = 3.86e26*10^(-0.4*(vj_ssp - 4.76)) / mass_stars_ssp(2999)
Lrc_ssp = 3.86e26*10^(-0.4*(rc_ssp - 4.76)) / mass_stars_ssp(2999)

$ b = burst strength in fraction of the mass of the main galaxy

num_iter = 4
time_line = fltarr(num_iter+1,5,2)

for i=0,num_iter do begin 

b = 0.008 + 0.003*i

$ Total luminosity, adding the underlying galaxy and the burst
Luj_tot = Luj_ug + b*Luj_ssp
Lbj_tot = Lbj_ug + b*Lbj_ssp
Lvj_tot = Lvj_ug + b*Lvj_ssp
Lrc_tot = Lrc_ug + b*Lrc_ssp

$ Calculate new magnitudes for the summed luminosity
$ Also, create a larger array with some pre-burst colors
Muj_tot = 4.76 - 2.5*alog10(Luj_tot/3.86e26)
Muj_tot = [uj_ug(300:step_ug), Muj_tot]
Mbj_tot = 4.76 - 2.5*alog10(Lbj_tot/3.86e26)
Mbj_tot = [bj_ug(300:step_ug), Mbj_tot]
Mvj_tot = 4.76 - 2.5*alog10(Lvj_tot/3.86e26)
Mvj_tot = [vj_ug(300:step_ug), Mvj_tot]
Mrc_tot = 4.76 - 2.5*alog10(Lrc_tot/3.86e26)
Mrc_tot = [rc_ug(300:step_ug), Mrc_tot]

$ Calculate new colors
UmB = Muj_tot-Mbj_tot
BmV = Mbj_tot-Mvj_tot
BmR = Mbj_tot-Mrc_tot
VmR = Mvj_tot-Mrc_tot
UmB_ssp = uj_ssp-bj_ssp
BmV_ssp = bj_ssp-vj_ssp
BmR_ssp = bj_ssp-rc_ssp
VmR_ssp = vj_ssp-rc_ssp

!x.thick=3
!y.thick=3
!p.charsize=1.5
!p.charthick=3.0
!p.thick=3.0

print, b

if (i EQ 0) then begin
	plot, BmV, UmB, xtitle='B-V', ytitle='U-B', xrange=[0.35,0.45], yrange=[-0.4,-0.05], xstyle=1, ystyle=1 


endif else begin
	oplot, BmV, UmB
endelse

time_line(i,0:4,0) = BmV(35:39)
time_line(i,0:4,1) = UmB(35:39)

; 0-33 are the underlying galaxy.  34 starts the burst colors, and 35 is 40Myr
times = ['40 Myr ', '80 Myr ', '110 Myr', '150 Myr', '180 Myr']
if b EQ 0.008 then for j=0,3 do label_symbol, BmV(j+35), UmB(j+35), 3, times(j), charsize=1.0

endfor

x_err = fltarr(14)
y_err = fltarr(14)
x_err(*) = 0.00
y_err(*) = 0.10
oploterror, BmV_ew(0:3), UmB_ew(0:3), x_err, y_err, psym=4
oploterror, BmV_ew(4:6), UmB_ew(4:6), x_err, y_err, psym=5
oploterror, BmV_ew(8:9), UmB_ew(8:9), x_err, y_err, psym=5
oploterror, BmV_ew(10:11), UmB_ew(10:11), x_err, y_err, psym=2
oploterror, BmV_ew(12), UmB_ew(12), x_err, y_err, psym=4
oploterror, BmV_ew(7), UmB_ew(7), x_err, y_err, psym=4
oplot, [0.7,0.723], [-0.7,-0.675]
for j=0,num_iter do oplot, time_line(*,j,0), time_line(*,j,1)
xyouts, 0.7, -0.75, textoidl('A_{\lambda}')
phot_regions = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N']

for i=0,3 do label_symbol, BmV_ew(i), UmB_ew(i), 4, phot_regions(i), charsize=1.0
for i=4,6 do label_symbol, BmV_ew(i), UmB_ew(i), 5, phot_regions(i), charsize=1.0
for i=8,9 do label_symbol, BmV_ew(i), UmB_ew(i), 5, phot_regions(i), charsize=1.0
label_symbol, BmV_ew(7), UmB_ew(7), 4, phot_regions(7), charsize=1.0
label_symbol, BmV_ew(10), UmB_ew(10), 2, phot_regions(10), charsize=1.0
label_symbol, BmV_ew(11), UmB_ew(11), 2, phot_regions(11), charsize=1.0
label_symbol, BmV_ew(12), UmB_ew(12), 4, phot_regions(12), charsize=1.0, /left

label_symbol, 0.15, 0.4, 4, ' Good U Magnitudes'
label_symbol, 0.15, 0.35, 2, ' Good U Magnitudes, but Strong Emission Lines'
label_symbol, 0.15, 0.3, 5, ' Low S/N in U'

device, /close
set_plot, 'x'


end

