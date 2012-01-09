
pro plot_focalplane, filename

    n = 17 ; the circle will be "created" with 17 data points (vertices)
    theta = findgen(n)/(n-1.0)*360.0*!DtoR ; 
    x = sin(theta)
    y = cos(theta)
    usersym, x, y

    x = read_dat(filename)
    x *= !radeg

    psfile_open, 'focalplane.eps', /thick, /encapsulated, /square

    plot, x[2,*], x[3,*], psym=8, symsize=0.5, $
        yrange=[-1.,0.1], xrange=[-1.1,-0.1], /xs, /ys, $
        xtitle='focalplane x [deg]', ytitle='focalplane y [deg]'

    n = n_elements(x[3,*])
    for i = 0,n-1 do begin
        x0 = x[2,i]
        y0 = x[3,i]
        dx = x[4,i]*0.0002
        dy = x[5,i]*0.0002
        oplot, [x0-dx,x0+dx], [y0-dy,y0+dy]
    endfor

    psfile_close

end

