; Code written by Elizabeth Wehner in IDL
;
; dE = dwarf Elliptical galaxy
; This project explores the relationship between the spheroidal component of dwarf
; elliptical galaxies with the mass of the system. Code below calculates the mass of the
; dE nuclei based on previously published M/L ratios. Results are then plotted alongside
; data for the mass of central black holes as compared to galaxy mass. Our results show
; a mass-continuum from the largest galaxies with supermassive black holes all the way down
; to tiny dwarf ellipticals, suggesting the spheroidal components of dEs are
; potentially thwarted black holes.
; 
; Results from the following code are published here:
; http://iopscience.iop.org/1538-4357/644/1/L17/
;


pro de

set_plot, 'ps'
device, filename='fa.eps', color=64, encapsul=1, /cmyk

!x.thick = 3.0
!y.thick = 3.0
!p.thick = 3.0
!p.charsize = 1.5
!p.charthick = 3.0
d = FSC_COLOR(/AllColors, ColorStructure=c)

$ Read in data from Lotz, Miller, & Fergus 2004 and Tremaine et al. 2002
readcol, 'lotz_table2.dat', format='A,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F', name, M_B, V, V_err, VmI, VmI_err, e, r0_arcsec, r0_err, r0_kpc, MV_n, MV_n_err, VmIn, VmIn_err, NGC, NGC_err
readcol, 'tremaine_table1.dat', format='A,A,F,F,F,F,F,F,F,F,A', t_name, type, t_M_B, M_BH, power, M_low, M_high, sigma1, D_mpc, MtoL, Band 

$ Read in new data (Gebhardt et al. 2003)
readcol, 'new_data.dat', format='A,A,F,F,F,F,F,F,F,F,F,F,F,F', new_name, new_type, new_M_B, new_M_BH, new_M_power, new_M_low, new_M_high, new_sigma, new_Dist, new_MtoL, new_a, new_b, new_c, new_d

for i=0,11 do print, new_name(i), new_M_B(i), new_M_BH(i)*10^(new_M_power(i))


M_BH = M_BH * 10^power
M_BH_err = (M_high-M_low) / 2.0 * 10^power / M_BH

$ Absolute V-magnitude of each dE's nucleus was used to calculate absolute L_V
$ which was then multiplied by a M/L ratio to compute the dEn_Mass
$ dEn_Mass = Mass (in solar units) of the nucleus in the dE
$ dE_L_B = B-band abs. luminosity of the entire dE, including the nucleus, in L_sun units
$ dE_L_Vn = V-band abs. magnitude of the nucleus only, in L_sun units

dE_L_B = 10^((5.48-M_B)/2.5)
dE_L_Vn = 10^((4.83-MV_n)/2.5)
;for i=0,37 do print, MV_n(i), dE_L_Vn(i), dE_L_V_err(i)

; t = tremaine, luminosity is for galaxy surrounding black hole.
t_L_B = 10^((5.48-t_M_B)/2.5)

MtoL_dEs = 3.0

M_BH_hosts = fltarr(30)
for j=0,29 do if (MtoL(j) NE -99.0) then M_BH_hosts(j) = t_L_B(j) * MtoL(j) else M_BH_hosts(j) = t_L_B(j) * 4.0

M_dE_nuc = dE_L_Vn*3.0

$ Mass to light ratio:

M_dE_gal = 3720.0*dE_L_B^(0.7)
dEn_Mass = dE_L_Vn * MtoL_dEs
dEn_Mass_err = (10^((4.83-(MV_n - MV_n_err))/2.5)*MtoL_dEs - dEn_Mass) / dEn_Mass


print, 'for the dE,N Galaxies: '
print, '(Luminosities in L_sun and Masses in M_sun)'
print, '   '
print, '        M_B     L_B          Mass_dE      Mass_Nuc '
for q=0,43 do print, M_B(q), dE_L_B(q), M_dE_gal(q), dEn_Mass(q)

print, 'for the BH Galaxies: '
print, '(Luminosities in L_sun and Masses in M_sun)'
print, '   '
print, '        M_B       L_B         Mass_host    Mass_BH '
for q=0,29 do print, M_B(q), t_L_B(q), M_BH_hosts(q), M_BH(q)

;for n=0,43 do print, dE_L_Vn(n), 10^((4.83-(MV_n(n) - MV_n_err(n)))/2.5), 10^((4.83-(MV_n(n) - MV_n_err(n)))/2.5) - dE_L_Vn(n) , dEn_Mass_err(n)

Mass_err_both = [dEn_Mass_err,M_BH_err]

dEn_weights = dEn_Mass / dEn_Mass
M_BH_weights = M_BH / M_BH


$ Plot Results: 

print, 'For all fits, Y-intercept is listed first, then slope' 
print, '  '

$ Black Holes only:
;plot, t_M_B, alog10(M_BH), psym=4, xtitle='M_B', ytitle='log (Mass_BH)', xrange=[-14,-22]
$ dE Nucleus mass versus abs B-luminosity
;plot, dE_L_B, alog10(dEn_Mass), psym=5, xtitle='L_B of dEs', ytitle='log(Mass_nucleus)
$ dE Nucleus mass (log) versus abs. B magnitude
;plot, M_B, alog10(dEn_Mass), psym=5, xtitle='M_B of total dE', ytitle='log(Mass_nucleus)


; ---- M_V of dE_nuc versus M_B dE Galaxy ----
; dE Nucleus Magnitude versus dE galaxy magnitude:
plot, M_B, MV_n, xtitle='M_B of dE Galaxy', ytitle='M_V of dE Nucleus', psym=4
result_grant = polyfitw(M_B, MV_n, dEn_weights, 1, yfit_grant, yband_grant, sigma_grant, corrm_grant)
;print, 'result_grant = ', result_grant
oplot, M_B, yfit_grant, lines=0
xyouts, -12.3, -12.5, 'Slope = 0.43'

device, /close

; ---- M_CMO vs. Luminosity of the Spheroid  ----

device, filename='f1.eps', color=64, encapsul=1, /cmyk
plot, alog10(t_L_B), alog10(M_BH), xtitle=textoidl('Log(L_{B,Sph})'), ytitle=textoidl('Log(M_{CMO})'), xrange=[6.0,11.0], yrange=[4.0,10.0], /nodata
oplot, alog10(t_L_B), alog10(M_BH), psym=4, color=c.red
oplot, alog10(dE_L_B), alog10(dEn_Mass) , psym=1, color=c.blue
legend,['dE Nuclei','Black Holes'],psym=[1,4], color=[c.blue,c.red]

t_L_B_weights = t_L_B / t_L_B
;print, t_L_B, M_BH
result_test = POLYFITW(alog10(t_L_B), alog10(M_BH), t_L_B_weights, 1, yfit_test, yband_test, sigma_test, corrm_test)
  ; print the fitted coefficients
  print, 'log(M_BH) versus log(t_L_B):  '
  print, result_test
  ; oplot the fitted function
  oplot, alog10(t_L_B(2:29)), yfit_test(2:29), lines=0
print,'laa ',  where(yfit_test GE 7.0), '  adsf  ', alog10(t_L_B), yfit_test

result_test2 = POLYFITW(alog10(dE_L_B), alog10(dEn_Mass), dEn_weights, 1, yfit_test2, yband_test2, sigma_test2, corrm_test2)
  ; print the fitted coefficients
  print, 'log(dEn_Mass) versus log(dE_L_B):'
  print, result_test2
  ; oplot the fitted function
  ;oplot, alog10(dE_L_B(12:43)), yfit_test2(12:43), lines=0
  ;print, where(yfit_test2 LE 6.5)

grant_slope_x=findgen(88)
grant_slope_x=grant_slope_x / 40.0 + 6.84
grant_slope_y=grant_slope_x * 0.7 + 0.73421
oplot, grant_slope_x, grant_slope_y, lines=0

device, /close

; ---- M_galaxy versus M_CMO ----



;m/l = 1
; for actual M_dE_gal, see above...
;M_dE_gal = 8110.0*dE_L_B^(0.7)
;M_dE_gal = 13300.0*dE_L_B^(0.7)
;M_dE_gal = 34.7*dE_L_B^(0.7) + dE_L_B*6.0
;M_dE_gal =  dE_L_B*6.0

;M_BH_err = ((M_low + M_high) / 2) * 10.0^power
;M_dE_nuc_err = dE_L_V_err * 3.0

device, filename='f2.eps', color=64, encapsul=1, /cmyk

plot, alog10(M_BH_hosts), alog10(M_BH), ytitle=textoidl('Log(M_{CMO})'), xtitle=textoidl('Log(M_{Sph})'), xrange=[8,12], yrange=[4,10], /nodata
oplot, alog10(M_BH_hosts), alog10(M_BH), psym=4, color=c.red
result_test3 = polyfitw(alog10(M_BH_hosts), alog10(M_BH), M_BH_weights, 1, yfit_test3, yband_test3, sigma_test3, corrm_test3)
   print, 'log(m_bh) versus log(m_cmo): ' 
   print, result_test3
   ;oplot, alog10(M_BH_hosts), yfit_test3, lines=0

legend,['dE Nuclei','Black Holes'],psym=[1,4], color=[c.blue,c.red]

oplot, alog10(M_dE_gal), alog10(M_de_nuc), psym=1, color=c.blue
result_test4 = polyfitw(alog10(M_dE_gal), alog10(M_de_nuc), dEn_weights, 1, yfit_test4, yband_test4, sigma_test4, corrm_test4)
   print, 'log(m_dEnuc) versus log(m_dE): '
   print, result_test4
   ;oplot, alog10(M_dE_gal), yfit_test4, lines=0

M_sph_both = [M_dE_gal, M_BH_hosts]
M_cmo_both = [M_dE_nuc, M_BH]
;M_err_both = [M_dE_nuc_err, M_BH_err]
M_err_both = M_cmo_both / M_cmo_both

;print, M_err_both
result_bothy = polyfitw(alog10(M_sph_both), alog10(M_cmo_both), M_err_both, 1, yfit_bothy, yband_bothy, sigma_bothy, corrm_bothy)
print, 'log(M_cmo) vs. log(M_sph) and sigma Y PolyfitW = '
print, result_bothy, 'sigma  = ', sigma_bothy, '  yband ', yband_bothy, ' corrm ', corrm_bothy
oplot, alog10(M_sph_both), yfit_bothy
plots, [8,9.7], [7,7], lines=1
plots, [9.7,9.7], [4,7], lines=1

temp_errs = M_err_both - 1
result_bothr = poly_fit(alog10(M_sph_both), alog10(M_cmo_both), 1, chisq=chisq_both, covar=covar_both, sigma=sigma_both, yfit=yfit_both)
print, 'log(M_cmo) vs. log(M_sph) and sigma Y = Polyfit'
print, result_bothr, 'sigma  = ', sigma_both, '  chisq=', chisq_both

result_bothx = poly_fit(alog10(M_cmo_both), alog10(M_sph_both), 1, chisq=chisq_bothx, covar=covar_bothx, sigma=sigma_bothx, yfit=yfit_bothx)
print, 'log(M_cmo) vs. log(M_sph) and sigma X = Polyfit'
print, result_bothx, 'sigma  = ', sigma_bothx,'  chisq= ', chisq_bothx

average_result = (result_bothx + result_bothr) / 2.0
;average_result = [average_result,(result_bothx(1) + result_bothy(0)) / 2.0]

yfitx_real = (yfit_bothx - result_bothx(0)) / result_bothx(1)
;oplot, result_bothx(1)*alog10(M_cmo_both)+result_bothx(0), yfitx_real, color=c.purple
;oplot, yfit_bothx, alog10(M_cmo_both), color=c.purple

x_avg = findgen(85)
x_avg = x_avg / 25.0 + 8.37
y_avg = 1.14*x_avg - 4.05 

oplot, x_avg, y_avg, lines=2

mx = result_bothx(1)
my1 = result_bothy(1)
bx = result_bothx(0)
by1 = result_bothy(0)

my2 = 1.0 / mx
by2 = -1*bx / mx

m_avg = (my2 + my1) / 2.0
b_avg = (by2 + by1) / 2.0

;print, alog10(M_cmo_both), yfitx_real, yfit_bothx
print, 'average result = "
print, 'm:', m_avg, '  b:', b_avg

device, /close
; ---- Mass of the CMO versus M_B ----


device, filename='fb.eps', color=64, encapsul=1, /cmyk
$ both Black Holes and dE Nuclei
plot, t_M_B, alog10(M_BH), psym=4, xtitle=textoidl('M_{B}'), ytitle=textoidl('log M_{CMO}'), xrange=[-11,-22], yrange=[4,10]
oplot, M_B, alog10(dEn_Mass), psym=1

$ --- fit line to dEns ---
result_dE = POLYFITW(M_B, alog10(dEn_Mass), dEn_weights, 1, yfit_dE, yband_dE, sigma_dE, corrm_dE)
  ; print the fitted coefficients
  print, 'log(dEn_Mass) vs. M_B: '
  print, result_dE
  ; oplot the fitted function
  oplot, M_B[6:43], yfit_dE[6:43], lines=0

;testing different fits: (ok - they work)
;result_dE = linfit(M_B, alog10(dEn_Mass),  yfit=yfit_dE)
;print, 'results from linfit: ', result_dE
;result_dE = poly_fit(M_B, alog10(dEn_Mass), 1, yfit=yfit_dE)
;print, 'results from polyfit: ', result_dE


$ --- fit line to BHs ---
result_BH = POLYFITW(t_M_B, alog10(M_BH), M_BH_weights, 1, yfit_BH, yband_BH, sigma_BH, corrm_BH)
  ; print the fitted coefficients
  print, 'log(M_BH) vs. t_M_B:  '
  print, result_BH
  ; oplot the fitted function
  oplot, t_M_B, yfit_BH, lines=0

x_intersect = (result_dE(0) - result_BH(0)) / (result_BH(1) - result_dE(1))
print, 'X intersection = ', x_intersect

legend,['dE Nuclei','Black Holes'],psym=[1,4], color=[c.blue,c.red]

;end


$ --- compare with alphas of slope from other papers ---

device, /close

set_plot, 'x'

end
