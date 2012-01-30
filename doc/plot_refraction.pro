
pro oplot_ref, x, temp, pressure, humidity, _extra=extra
    iwh = where( abs(x[0,*]-temp) lt 1e-6 and abs(x[1,*]-pressure) lt 1e-6 and abs(x[2,*]-humidity) lt 1e-6 )
    oplot, x[4,iwh], x[5,iwh], _extra=extra
end

pro plot_refraction, filename

    x = read_dat(filename, nskip=1)

    psfile_open, 'refraction.eps', /encapsulated, /golden, /thick
    color = getcolor()
    plot, [35,60], [20,70], /nodata, $
        xtitle='altitude [deg]', ytitle='refraction [arcsec]'

    oplot_ref, x, 20., 548., 1.0, color=color.orange
    oplot_ref, x, 20., 560., 1.0, color=color.orange, linestyle=1

    oplot_ref, x, -20., 548., 1.0, color=color.green
    oplot_ref, x, -20., 560., 1.0, color=color.green, linestyle=1

    oplot_ref, x, -20., 548., 0.0, color=color.blue
    oplot_ref, x, -20., 560., 0.0, color=color.blue, linestyle=1

    oplot_ref, x, 20., 548., 0.0, color=color.red
    oplot_ref, x, 20., 560., 0.0, color=color.red, linestyle=1

    legend, ['"hot & humid" (20 C, 100%)', '"cold & humid" (-20 C, 100%)', $
        '"cold & dry" (-20 C, 0%)', '"hot & dry" (20 C, 0%)'], $
        color=[color.orange, color.green, color.blue, color.red], $
        linestyle=[0,0,0,0], charsize=1.0, /bottom, /left, box=0
    legend, ['low pressure (548 mbar)', 'high pressure (560 mbar)'], linestyle=[0,1], $
        /top, /right, box=0, charsize=1.0

    psfile_close
end

