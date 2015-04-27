; Written by Elizabeth Wehner in IDL
; analysis of globular star cluster system in the nearby(ish) 
; cD massive elliptical galaxy, NGC 3311, and close companion galaxy
; of NGC 3309.
;
; Note! Code was written and adapted as needed over year-long analysis
; of intergalactic system. Therefore, code is not optimized for efficiency,
; but rather its creation was guided by the flow of the astronomical research.
;
; Public article available online: http://www.gemini.edu/node/247
; Professional journal article here:  http://iopscience.iop.org/0004-637X/681/2/1233/

pro mp

set_plot, 'ps'
!x.thick = 4.0
!y.thick = 4.0
!p.thick = 4.0
!p.charsize = 1.5
!p.charthick = 4.0
!p.multi=0

d = FSC_COLOR(/AllColors, ColorStructure=c)

;readcol, 'final_star_list.dat', x, y, i, i_err, chi_i, g, g_err, chi_g, dr
;readcol, 'calib.dat', x, y, i, i_err, chi_i, g, g_err, chi_g, gmi, gmi_err
readcol, 'cull.dat', x, y, i, i_err, chi_i, g, g_err, chi_g, gmi, gmi_err, jnk

title_irange = ''

; moved above plot of gmi vs i in order to examine limits for kmm
; SHOULD PROBABLY MOVE BACK TO AFTER PLOT gmi vs. i FOR FINAL PLOT
; correct for reddening:  E(g-i) = 0.158
e_gmi = 0.158
gmi = gmi-e_gmi

; Figure out which of those GCs are outside the noisy radii in each galaxy and 
; which are on the side of NGC 3311 -- CALCULATE good_gcs and one_sided

x_centre = 1037.0
y_centre = 1168.0
companion_x_centre = 796.0
companion_y_centre = 521.0
good_gcs = 0 
one_sided = 0

sz = size(i)
s = sz(1)
for l=0,s-1 do begin
  x_line = (y(l) - 1554.9)/(-0.3731)
  ; also had (x(l) GE x_line) in if statement below... (re-added 4.30.07)
  d_from_center = sqrt((x(l)-x_centre)^2 + (y(l)-y_centre)^2)
  d_from_companion = sqrt((x(l)-companion_x_centre)^2 + (y(l)-companion_y_centre)^2)
  if ((d_from_center GE 120.0) AND (d_from_companion GE 120.0)) then good_gcs=[good_gcs,l]
  if ((d_from_center GE 120.0) AND (d_from_companion GE 120.0) AND (x(l) GE x_line)) then one_sided=[one_sided,l]
end

; change over i, gmi to only those on n3311 side of chip...
; good_gcs represents all clusters except those in an inner circle
; around each of NGC 3311 and NGC 3309.  One_sided excludes those 
; circles as well as all GCs on the side of NGC 3309.  One_sided is
; to be used to find the radial distribution of the clusters, in 
; order to avoid contamination from the 3309 clusters.

sz = size(good_gcs)
good_gcs = good_gcs(1:sz(1)-1)
;print, 'testing   ', good_gcs
i = i(good_gcs)
g = g(good_gcs)
x = x(good_gcs)
y = y(good_gcs)
gmi = gmi(good_gcs)

result23 = where(i LE 23.0)
result24 = where(i LE 24.0)
num24 = size(result24)
num24 = num24(1)
num23 = size(result23)
num23 = num23(1)
print, 'Number greater than 23: ', num23, '  24:  ', num24
nblends23 = num23^2/2.0 * (3.14159*0.5^2/(5.5*60)^2)
nblends24 = num24^2/2.0 * (3.14159*0.5^2/(5.5*60)^2)

print, 'nblends 23, 24:  ', nblends23, nblends24

device, /inches, xsize=6, filename='xyplot.eps'
plot, x, y, xtitle='X (pixels)', ytitle='Y (pixels)', psym=3

; output list of colors

sz = size(i)
s = sz(1)
free_lun, 1
openw, 1, 'gmi_list.dat'
for l=0,s-1 do if ((i(l) le 23.0) AND (i(l) ge 22.0) AND (gmi(l) LE 1.3) AND (gmi(l) GE 0.20)) then printf, 1, gmi(l)
close, 1
free_lun, 1

plotsym, 0, /fill
; Plot g-i vs. i

bright_lim = 22.15
blue_lim = 0.60
red_lim = 1.1

device, /inches, color=64, xsize=7, ysize=7, encapsulated=1, filename='gmi.eps'
plot, gmi, i, xrange=[-1,2], yrange=[28,20.5], psym=8, symsize=0.3, xstyle=1, ystyle=1, xtitle=textoidl("g\'-i\'"), ytitle = textoidl("i\'")
plots, [blue_lim,blue_lim],[bright_lim,20.5], color=c.red, thick=2
plots, [blue_lim,1.1],[20.5,20.5], color=c.red, thick=2
plots, [1.1,1.1],[bright_lim,20.5], color=c.red, thick=2
plots, [blue_lim,1.1],[bright_lim,bright_lim], color=c.red, thick=2
device, /close


; convert x,y coordinates of globular clusters into their linear coord along n3311-n3309 line

new_coord_subset = where(i le 24.5 AND gmi LE 1.1 AND gmi GE 0.4)
y_ncsub = y(new_coord_subset)
x_ncsub = x(new_coord_subset)
r_ncsub = sqrt((x_ncsub-x_centre)^2 + (y_ncsub-y_centre)^2)
sz = size(x_ncsub)
s = sz(1)
line_coord = fltarr(s)
; b = y intercept, par = parallel line, per = perpendicular
b_par = -1612.3
m_par = 2.68
m_per = -0.3731
for h=0,s-1 do begin
    b_per = y_ncsub(h) - m_per*x_ncsub(h)
    x_line = (b_par - b_per) / (m_per - m_par)
    y_line = m_per*x_line + b_per
    x_yeqzero = (y_line - b_par) / m_par
    line_coord(h) = sqrt((x_line-x_yeqzero)^2 + y_line^2)
end

line_coord_hist = histogram(line_coord, binsize=100.0, locations = line_coord_hist_locations)
device, encapsulated=1, filename='new_coords.eps'
plot, line_coord_hist_locations, line_coord_hist, psym=4, ytitle='Number', xtitle='Distance'
device, /close

free_lun, 4
openw, 4, 'new_coords.dat'
sz = size(line_coord)
printf, 4, '     coord '
for j=0,sz(1)-1 do printf, 4, line_coord(j)
close, 4

free_lun, 5
openw, 5, 'testcoords.dat'
sz = size(line_coord_hist)
for j=0,sz(1)-1 do printf, 5, line_coord_hist_locations(j), line_coord_hist(j)
close, 5

; isolate brightest clusters and plot their x, y locations

; SUBSELECT THOSE IN RED BOX IN GMI.EPS PLOT??????

bright_lim = 22.15
blue_lim = 0.60
red_lim = 1.1

bright = where(i LE bright_lim AND gmi LE red_lim AND gmi GE blue_lim, count1)
; print, 'count1 = ', count1
i_bright=i(bright)
g_bright=g(bright)
gmi_bright = gmi(bright)
x_bright=x(bright)
y_bright=y(bright)
; print, i_bright(sort(i_bright))

free_lun, 3
openw, 3, 'bright_coords.dat'
printf, 3, '# Region file format: DS9 version 4.0'
printf, 3, '# Filename: /Users/wehnere/final_dec06/i.fits'
printf, 3, 'global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source'
printf, 3, 'image'
for q=0,count1-1 do printf, 3, 'circle(', trim(x_bright(q),2), ',', trim(y_bright(q),2), ',12)'
close, 3
device, /inches, xsize=6, filename='bright_objs.eps'
plot, x(bright), y(bright), xtitle='X (pixels)', ytitle='Y (pixels)', psym=2

device, /inches, xsize=6, filename='xyplot_orig.eps'
plot, x, y, xtitle='X (pixels)', ytitle='Y (pixels)', psym=3




; CALCULATE THE RADIAL DISTRIBUTION USING one_sided

sz = size(one_sided)
one_sided = one_sided(1:sz(1)-1)
i_one_sided = i(one_sided)
x_one_sided = x(one_sided)
y_one_sided = y(one_sided)
gmi_one_sided = gmi(one_sided)

my_bin_size = 300.0
num_bins = 5.0
gc_errors = fltarr(num_bins)
dgto_errors = fltarr(num_bins)
sqkpc_per_pix = 0.001439

r_one_sided = sqrt((x_one_sided-x_centre)^2 + (y_one_sided-y_centre)^2)
one_sided_hist = histogram(r_one_sided, binsize=my_bin_size, nbins=num_bins, locations = one_sided_hist_locations)
area_circles = fltarr(num_bins)
area_rings = fltarr(num_bins)

print, 'one sided hist:', one_sided_hist(0), one_sided_hist(1)
help, one_sided_hist_locations

; THIS IS THE CORRECT AREA CODE:  (the following loop)
; since i'm done making the radial distribution plot, i'm replacing it with the old
; area code (which follows) because it's much less time intensive.  Restore this code 
; and comment out the following 3 lines if you want the correct radial distribution plot
; (modified May 31, 2007)

for q=0,num_bins-1 do begin
	r1 = q*my_bin_size
	r2 = (q+1)*my_bin_size
	for x_pix = 1,2261 do begin
		for y_pix = 1,2260 do begin
			r_pix = sqrt((x_pix-x_centre)^2 + (y_pix - y_centre)^2)
			if ((r_pix GE r1) AND (r_pix LE r2) AND (r_pix GE 120.0)) then area_rings(q) = area_rings(q) + 1
		end
	end
end

; the following three lines are NOT really the correct area calculations
; but are much faster.  leave this in for speed, take them out (and replace w/above)
; for the radial distribution of the DGTOs and GCs plot
;
;for j=1,num_bins do area_circles(j-1) = 3.14*(j*my_bin_size)^2
;area_rings(0)= area_circles(0)
;for j=1,num_bins-1 do area_rings(j) = area_circles(j)-area_circles(j-1)

area_rings = area_rings * sqkpc_per_pix

r_bright = sqrt((x_bright-x_centre)^2 + (y_bright-y_centre)^2)
device, /inches, xsize=7, ysize=5, filename='bright_rad_distrib.eps'
bright_hist = histogram(r_bright, binsize=my_bin_size, nbins=num_bins, locations = bright_hist_locations)
dgto_errors = sqrt(bright_hist)/area_rings
help, bright_hist

; bright_fit = gaussfit(bright_hist_locations, bright_hist/area_rings, my_coeffs)
; print, 'my coeffs:  ', my_coeffs

bright_fit = poly_fit(bright_hist_locations, bright_hist/area_rings, 4, yfit=bright_fit_points)
one_sided_fit = poly_fit(one_sided_hist_locations, one_sided_hist/area_rings, 2, yfit=one_sided_fit_points)
kpcperpix = 0.03849

plot, bright_hist_locations*kpcperpix, bright_hist/area_rings, psym=4, yrange=[0,bright_hist(0)/area_rings(0)], ytitle='Number per Square Kpc', xtitle='Radius (in Kpc)'
oplot, bright_hist_locations*kpcperpix, bright_fit_points, linestyle=1
norm_factor = one_sided_hist(0) / bright_hist(0) 
;norm_factor = 1.0/(total(bright_hist)/ total(one_sided_hist))
print, 'norm_factor:  ', norm_factor
gc_errors = sqrt(one_sided_hist)/area_rings / norm_factor
; print, gc_errors, dgto_errors, norm_factor
oplot, one_sided_hist_locations*kpcperpix, one_sided_hist/area_rings/norm_factor, psym=1
oplot, one_sided_hist_locations*kpcperpix, one_sided_fit_points/norm_factor, linestyle=3
legend, ['UCD Fit','GCS Fit (Normalized)'], linestyle=[1,3], /right
legend, ['UCDs','All GCs (Normalized)'],psym=[4,1], pos=[58.5,0.032], /right
oploterr, one_sided_hist_locations*kpcperpix, one_sided_hist/area_rings/norm_factor, gc_errors, 3
oploterr, bright_hist_locations*kpcperpix, bright_hist/area_rings, dgto_errors, 1


; print, one_sided_hist(0)/area_rings(0), norm_factor
; print, dgto_errors, area_rings
device, /close

num_ucds = size(r_bright, /dimensions)
cm_ucd_dist = fltarr(num_ucds)
cm_num_ucds = bright_hist(0)
cm_ring_area = area_rings(0)
cm_ucd_dist(0) = cm_num_ucds / total(bright_hist)

;for q=1,num_bins-1 do begin
;	cm_num_ucds += bright_hist(q)
;	cm_ring_area += area_rings(q)
;	cm_ucd_dist(q) = float(cm_num_ucds) / num_ucds
;	print, cm_num_ucds, num_ucds, cm_ucd_dist(q)
;end

result = sort(r_bright)
r_bright_temp = r_bright(result)
cm_ucd_dist(0) = 1.0
for q=1,28 do cm_ucd_dist(q) = cm_ucd_dist(q-1) + 1.0
cm_ucd_dist = cm_ucd_dist / 29.0
print, cm_ucd_dist

result = sort(r_ncsub)
r_ncsub_temp = r_ncsub(result)
num_gcs = size(r_ncsub_temp, /dimensions)
print, 'num_gcs =  ', num_gcs
num_one_sided = size(x_one_sided, /dimensions)
cm_gc_dist = fltarr(num_gcs)
cm_num_gcs = one_sided_hist(0)
cm_ring_area = area_rings(0)
cm_gc_dist(0) = cm_num_ucds / total(one_sided_hist)

num_gcs = fix(num_gcs(0))
result = sort(r_ncsub)
r_ncsub_temp = r_ncsub(result)
cm_gc_dist(0) = 1.0
num_gcs = size(r_ncsub_temp, /dimensions)
num_gcs = fix(num_gcs(0))
for q=1,num_gcs-1 do cm_gc_dist(q) = cm_gc_dist(q-1) + 1.0
cm_gc_dist = cm_gc_dist / num_gcs
print, cm_gc_dist

;for q=1,num_bins-1 do begin
;	cm_num_gcs += one_sided_hist(q)
;	cm_ring_area += area_rings(q)
;	cm_gc_dist(q) = cm_num_gcs / total(one_sided_hist)
;	print, cm_num_gcs, total(one_sided_hist), cm_gc_dist(q)
;end

device, filename='test.ps'
plot, r_bright_temp*kpcperpix, cm_ucd_dist, linestyle=1, ytitle=textoidl("Fraction within Radius R"), xtitle="Radius (in Kpc)"
oplot, r_bright_temp*kpcperpix, cm_ucd_dist, psym=4
oplot, r_ncsub_temp*kpcperpix, cm_gc_dist, linestyle=2
legend, ['UCDs','GCs'],linestyle=[1,2], pos=[57,0.2], /right
device, /close





; CALCULATE THE MASSES AND LUMINOSITIES OF THE DGTOs

dist_mod = 33.68

sz = size(r_bright)
count = 0
i_temp = fltarr(sz(1))
g_temp = fltarr(sz(1))
r_bright_temp = fltarr(sz(1))
x_bright_temp = fltarr(sz(1))
y_bright_temp = fltarr(sz(1))
gmi_temp = fltarr(sz(1))
M_g = fltarr(sz(1))
L_g = fltarr(sz(1))
Mass_bright = fltarr(sz(1))
for q=0,sz(1)-1 do if (gmi_bright(q) GE blue_lim) AND (gmi_bright(q) LE 1.1) then begin
 	M_g(count) = g_bright(q) - dist_mod
	L_g(count) = 10^((5.06-M_g(count))/2.5) 
	Mass_bright(count) = L_g(count) * 3.0
	i_temp(count) = i_bright(q)
	g_temp(count) = g_bright(q)
	r_bright_temp(count) = r_bright(q)
	x_bright_temp(count) = x_bright(q)
	y_bright_temp(count) = y_bright(q)
	gmi_temp(count) = gmi_bright(q)
	count = count + 1
	; print, count, '    ', trim(M_g,1), '   ', L_g, Mass_bright, i_bright(q), gmi_bright(q)
end

result = sort(i_temp)
i_temp = i_temp(result)
g_temp = g_temp(result)
gmi_temp = gmi_temp(result)
x_bright_temp = x_bright_temp(result)
y_bright_temp = y_bright_temp(result)
r_bright_temp = r_bright_temp(result)
kpc_per_pix = 0.03849
r_bright_temp = r_bright_temp*kpc_per_pix
M_g = M_g(result)
L_g = L_g(result)
Mass_bright = Mass_bright(result)
sz = size(i_temp)
print, '  '

 for q=0,sz(1)-1 do print, q+1, x_bright_temp(q), y_bright_temp(q), '    ', trim(M_g(q),1), '   ', trim(Mass_bright(q)/1e6,1), '     ',  trim(i_temp(q),2), '    ', trim(gmi_temp(q),2), '   ', trim(r_bright_temp(q),1)

; Separate out GCs by magnitude bin, and then explore their color properties

color_step = 0.1
metal_step = 0.2
magstep = 0.5
min_mag = 22
max_mag = 26.5

num_steps = (max_mag - min_mag) / magstep
avg_color = fltarr(num_steps)
red_avg_color = fltarr(num_steps)
blue_avg_color = fltarr(num_steps)
avg_metal = fltarr(num_steps)
middle_mag = fltarr(num_steps)

!p.multi=[0,3,2]
device, /inches, xsize=10, filename='color_bins.eps'
for k = 0,(num_steps-1) do begin
	lower_lim =  min_mag+k*magstep
	upper_lim = min_mag+(k+1)*magstep
	middle_mag(k) = (lower_lim + upper_lim) / 2.0
	mag_subset = where(i GT lower_lim AND i LE upper_lim, count)
	if count GT 0 then begin
		color_hist = histogram(gmi(mag_subset), binsize=color_step, locations = color_bin_locations)
	 	gc_pop = where(gmi(mag_subset) GT 0.4 AND gmi(mag_subset) LT 1.1)
		avg_color(k) = avg(gmi(mag_subset(gc_pop)))
		color_label = textoidl("<g-i>_0") + ' = ' + trim(avg_color(k),2)
		title_irange = trim(lower_lim,1) + ' < i < ' + trim(upper_lim,1)
		plot, color_bin_locations, color_hist, title = title_irange, psym=10, xrange=[-1,2], ytitle='Number', xtitle='g-i', yrange=[0,max(color_hist)*1.2], ystyle=1
		sz = size(color_hist)
		s = sz(1)
		plots, [color_bin_locations(s-1),color_bin_locations(s-1)],[0,color_hist(s-1)]
		plots, [color_bin_locations(0),color_bin_locations(0)],[0,color_hist(0)]
		max_color =  max(color_hist)
		xyouts, -0.9, max_color*1.1, color_label, charsize=0.8

		; find avg of red half, and avg of blue half;
		
		red_half = where(gmi(mag_subset) GT 0.8 AND gmi(mag_subset) LT 1.2)
		red_avg_color(k) = avg(gmi(mag_subset(red_half)))
		red_fit = poly_fit(middle_mag, red_avg_color, 1, yfit=red_fit_points)
		blue_half = where(gmi(mag_subset) GT 0.3 AND gmi(mag_subset) LE 0.8)
		blue_avg_color(k) = avg(gmi(mag_subset(blue_half)))
		blue_fit = poly_fit(middle_mag, blue_avg_color, 1, yfit=blue_fit_points)
		;print, middle_mag(k), red_avg_color(k), blue_avg_color(k)
	endif
end
device, /close

!p.multi=[0,1,1]
device, /inches, xsize=7, ysize=5, filename='color_trend.eps'
plot, middle_mag, red_avg_color, xrange=[27,22], ytitle=textoidl("<g'-i'>_0"), xtitle=textoidl("i'_0"), psym=4, title='NGC 3311'
oplot, middle_mag, blue_avg_color, psym=5
oplot, middle_mag, blue_fit_points
oplot, middle_mag, red_fit_points
device, /close

!p.multi=[0,3,2]

; most of this metallicity stuff is crap -- we decided that the zero points of the 
; metallicity models in g' and i' were not well enough known (waiting on Mike West
; and his study of Milky Way GCs in g' and i')

metallicity = -7.531 + 9.577*gmi - 2.779*gmi^2
gc_low = -4.0
low_metal = -4.0
gc_high = 1.0
high_metal = 1.0

;print, low_metal, high_metal, gc_low, gc_high

device, /inches, xsize=10, filename='metal_bins.eps'
for k = 0,(num_steps-1) do begin
        lower_lim =  min_mag+k*magstep
        upper_lim = min_mag+(k+1)*magstep
        mag_subset = where(i GE lower_lim AND i LE upper_lim)
        if mag_subset(0) GT 0 then begin
                metal_hist = histogram(metallicity(mag_subset), binsize=metal_step, locations = metal_bin_locations)
		gc_pop = where(metallicity(mag_subset) GT gc_low)
		avg_metal(k) = avg(metallicity(mag_subset(gc_pop)))
		metal_label = '<m/H> = ' + trim(avg_metal(k),2)
                title_irange = trim(lower_lim,1) + ' < i < ' + trim(upper_lim,1)
                plot, metal_bin_locations, metal_hist, title = title_irange, psym=10, ytitle='Number', xtitle='m/H', yrange=[0,max(metal_hist)*1.2], xrange=[low_metal,high_metal], ystyle=1
                sz = size(metal_hist)
                s = sz(1)
                plots, [metal_bin_locations(s-1),metal_bin_locations(s-1)],[0,metal_hist(s-1)]
                plots, [metal_bin_locations(0),metal_bin_locations(0)],[0,metal_hist(0)]
		max_metal =  max(metal_hist)
		xyouts, -3.8, max_metal*1.1, metal_label, charsize=0.8
        endif
end
device, /close

!p.multi=[0,1,1]
device, /inches, xsize=7, filename='avg_color.eps'
above_lim = where(middle_mag LE 25.0)
;print, middle_mag(above_lim)
plot, middle_mag, avg_color, ytitle=textoidl("<g-i>_0"), xtitle=textoidl("i\'"), psym=4, xrange=[22,26]
plots, [25.0,25.0],[0.0, 1.5], linestyle=1
sig = stddev(avg_color(above_lim))
fit = linfit(middle_mag(above_lim), avg_color(above_lim), yfit = yfit_color)
;oplot, middle_mag(above_lim), yfit_color
device, /close

device, filename='avg_metal.eps'
plot, middle_mag, avg_metal, ytitle='<m/H>', xtitle='i', psym=4, yrange=[-4.0,0.0], xrange=[22,26]
fit = linfit(middle_mag(above_lim), avg_metal(above_lim), yfit=yfit_metal)
plots, [25.0,25.0],[-4.0,0.0], linestyle=1
;oplot, middle_mag(above_lim), yfit_metal
device, /close

set_plot, 'x'

end
