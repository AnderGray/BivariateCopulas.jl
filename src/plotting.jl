######
# This file is part of the BivariateCopulas.jl package.
#
#   Plotting joint cdfs, densities and samples
#
#           University of Liverpool, Institute for Risk and Uncertainty
#
#                                           Author: Ander Gray
#                                           Email:  ander.gray@liverpool.ac.uk
######


wait_for_key() = (print(stdout, "press enter to continue"); read(stdin, 1); nothing);

function samplePlot(C :: AbstractCopula, N = 1)     # Generates a plot showing conditional sampling

    x = rand(N);    y = rand(N);
    ux = x;         uy = zeros(N);
    m = n;
    e = 1/m;        ygrid = range(0,stop=1,length = m);

    if !ismissing(C.func) e = sqrt(eps()); end

    for i=1:N
        conditional = (C(x[i] + e/2,ygrid) - C(x[i] - e/2,ygrid))./e
        uy[i] = findlast(conditional .<= y[i]) / m;
        println("-------------------------------------------")
        println("N: $i.")
        println("Independent ~U(0,1)'s:")
        println(hcat(x[i],y[i]))
        println();

        plot2(C.cdf,conditional,ux[i], uy[i], y[i])

        println("Samples from copula:");
        println(hcat(ux[i],uy[i]))
        println("-------------------------------------------")

        wait_for_key()
    end
    return hcat(ux,uy);
end

plotDensity(x :: AbstractCopula) = plot1(x.density);
plotCdf(x::AbstractCopula) = plot1(x.cdf);

function plotCdf(j::AbstractJoint; pn = 50) 

    A = j.cdf
    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    nm = round(m/ppn);
    A = A[1:Int(nm):end,1:Int(nm):end]
    plot3(A, j.xRange[1:Int(nm):end],j.yRange[1:Int(nm):end])

end

function plotDensity(j::AbstractJoint) 

    A = j.density
    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    nm = round(m/ppn);
    A = A[1:Int(nm):end,1:Int(nm):end]
    plot3(A, j.xRange[1:Int(nm):end],j.yRange[1:Int(nm):end])

end

plot(j::AbstractJoint) = plot4(j)
plotDen(j::AbstractJoint) = plot4(j,false)

plotContourCdf(j::AbstractJoint) = plot5(j)
plotContourDen(j::AbstractJoint) = plot5(j,false)




function plot(x :: AbstractCopula;title = "SurfacePlots", pn = 50, fontsize=18)
    A = x.cdf
    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0,stop = 1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    z = A[1:Int(nm):end,1:Int(nm):end]

    fig = figure(title,figsize=(10,10))
    ax = fig.add_subplot(1,1,1,projection="3d")
    #ax = fig.add_subplot(2,1,1)
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, alpha=1, linewidth=2, cmap=ColorMap("RdGy"))
    ax.view_init(45-27, 180+ 26)
    xlabel("y",fontsize=fontsize)
    ylabel("x",fontsize=fontsize)
    ax.zaxis.set_rotate_label(false);
    zlabel("C(x,y)", rotation = 0,fontsize=fontsize)
    #xticks(fontsize = fontsize); yticks(fontsize = fontsize);
    #PyPlot.title(title,fontsize=fontsize)
    tight_layout()

end

function plotDen(x :: AbstractCopula; title = "SurfacePlots", pn = 50, fontsize=18)
    A = x.density
    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0,stop = 1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    z = A[1:Int(nm):end,1:Int(nm):end]

    fig = figure(title,figsize=(10,10))
    ax = fig.add_subplot(1,1,1,projection="3d")
    #ax = fig.add_subplot(2,1,1)
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, alpha=1, linewidth=2, cmap=ColorMap("RdGy"))
    ax.view_init(45-27, 180+ 26)
    xlabel("y",fontsize=fontsize)
    ylabel("x",fontsize=fontsize)
    ax.zaxis.set_rotate_label(false);
    zlabel("C(x,y)", rotation = 0,fontsize=fontsize)
    #xticks(fontsize = fontsize); yticks(fontsize = fontsize);
    #PyPlot.title(title,fontsize=fontsize)
    tight_layout()

end

function plotContourDen(x; title = "SurfacePlots", pn = pn, fontsize=18)
    A = x.density
    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0,stop = 1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    z = A[1:Int(nm):end,1:Int(nm):end]
    fig = figure(title,figsize=(10,10))

    cp = contour(xgrid, ygrid, z,cmap=ColorMap("coolwarm"),levels = 20)
    xlabel("y",fontsize=fontsize)
    ylabel("x",fontsize=fontsize)

end

function plotContourCdf(x :: AbstractCopula; title = "SurfacePlots",fontsize=18)
    A = x.cdf
    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0,stop = 1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    z = A[1:Int(nm):end,1:Int(nm):end]
    fig = figure(title,figsize=(10,10))

    cp = contour(xgrid, ygrid, z,cmap=ColorMap("coolwarm"),levels = 20)
    xlabel("y",fontsize=fontsize)
    ylabel("x",fontsize=fontsize)

end



function plot1(A :: Array{<:Real,2})

#   https://github.com/gizmaa/Julia_Examples/blob/master/pyplot_surfaceplot.jl
    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0,stop=1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = round(m/ppn);

    z = A[1:Int(nm):end,1:Int(nm):end]

    fig = figure("SurfacPlots",figsize=(10,10))
    ax = fig.add_subplot(2,1,1,projection="3d")
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, alpha=0.8, linewidth=0.25,cmap=ColorMap("coolwarm"))
    xlabel("X")
    ylabel("Y")
    ax.set_zlim3d(0, max(maximum(A),1)) 
    PyPlot.title("Surface Plot")

    subplot(212)
    ax = fig.add_subplot(2,1,2)
    cp = contour(xgrid, ygrid, z,cmap=ColorMap("coolwarm"),levels = 15)
    ax.clabel(cp, inline=1, fontsize=10)
    xlabel("X")
    ylabel("Y")
    #ax.set_zlim3d(0, max(maximum(A),1))
    PyPlot.title("Contour Plot")
    tight_layout()

end

function plot2(A :: Array{<:Real,2}, cond:: Array{<:Real,1}, ux = missing, uy = missing, uuy = missing)

    if !isempty(get_fignums()) PyPlot.clf(); end

    m = size(A)[1];
    if m < pn; ppn = m; else ppn = pn; end

    x = y = range(0,stop=1,length=ppn)
    xgrid = repeat(x',ppn,1)
    ygrid = repeat(y,1,ppn)

    nm = m/ppn;

    z = A[1:Int(nm):end, 1:Int(nm):end]
    conditional = cond[1:Int(nm):end];

    fig = figure("ConditionalPlots",figsize=(10,10))
    ax = fig.add_subplot(2,1,1)
    cp = contour(xgrid, ygrid, z)
    ax.clabel(cp, inline=1, fontsize=10)
    xlabel("X")
    ylabel("Y")

    if !ismissing(ux)
        xx = ones(pn) * ux
        PyPlot.plot(xx,x, color = "red",linewidth=2,label = "Rand1: $ux");
    end

    subplot(212)
    ax = fig.add_subplot(2,1,2)
    PyPlot.plot(y, conditional, color = "black")
    if !ismissing(uy) && !ismissing(uuy)
        nums1 = Int(round(pn*uy))
        nums2 = Int(round(pn*uuy))
        yy = ones(nums1) * uuy
        uyy = ones(nums2) * uy
        PyPlot.plot(y[1:nums1],yy,color = "red",linewidth=2, label= "Rand2: $uy");
        PyPlot.plot(uyy,y[1:nums2],color = "blue",linewidth=2);
    end
    xlabel("Y")
    ylabel("conditional CDF")
    tight_layout()
end

function plot3(A :: Array{<:Real,2}, x = range(0,1,length = pn), y = range(0,1,length = pn))    # Won't work if n!=200

#   https://github.com/gizmaa/Julia_Examples/blob/master/pyplot_surfaceplot.jl

    #m1 = size(A)[1];
    #m2 = size(A)[2];
    #if m1 < pn; ppn1 = m1; else ppn1 = pn; end      # ppn1 is the number of points to plot in X
    #if m2 < pn; ppn2 = m2; else ppn2 = pn; end      # ppn2 is the number of points to plot in Y

    #nm = round(m/ppn);                              # The descritization of the domain

    #x = range(rangeX[1],rangeX[2],length=pn)
    #y = range(rangeY[1],rangeY[1],length=pn)
    xgrid = repeat(x',pn,1)
    ygrid = repeat(y,1,pn)
    z = A'

    fig = figure("SurfacePlots",figsize=(10,10))
    ax = fig.add_subplot(2,1,1,projection="3d")
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("coolwarm"), alpha=0.8, linewidth=0.25)
    #plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("Gray"), alpha=0.8, linewidth=0.25)
    xlabel("X")
    ylabel("Y")
    PyPlot.title("Surface Plot")

    #subplot(212)
    ax = fig.add_subplot(2,1,2)
    cp = contour(xgrid, ygrid, z, cmap=ColorMap("coolwarm"))
    ax.clabel(cp, inline=1, fontsize=10)
    xlabel("X")
    ylabel("Y")
    PyPlot.title("Contour Plot")
    tight_layout()

end

function plot4(J :: AbstractJoint, CDF = true)    # Won't work if n!=200

#   https://github.com/gizmaa/Julia_Examples/blob/master/pyplot_surfaceplot.jl

    x = J.xRange;       y = J.yRange;
    M1 = J.marginal1;   M2 = J.marginal2;

    if CDF; z = J.cdf; title = "cdf" else z = J.density; title = "pdf";end

    m = size(z)[1];
    if m < pn; ppn = m; else ppn = pn; end

    nm = round(m/ppn);
    z = z[1:Int(nm):end,1:Int(nm):end]

    x = x[1:Int(nm):end]
    y = y[1:Int(nm):end]

    xgrid = repeat(x',pn,1)
    ygrid = repeat(y,1,pn)


    fig1 = figure("SurfacePlots",figsize=(10,10))
    #ax = fig.add_subplot(2,1,1,projection="3d")
    plot_surface(xgrid, ygrid, z, rstride=2,edgecolors="k", cstride=2, cmap=ColorMap("coolwarm"), alpha=0.8, linewidth=0.25)
    xlabel("X")
    ylabel("Y")

    PyPlot.title("$title surface plot")

    tight_layout()

end

function plot5(J :: AbstractJoint, CDF = true)    # Won't work if n!=200

    #   https://github.com/gizmaa/Julia_Examples/blob/master/pyplot_surfaceplot.jl
    
        x = J.xRange;       y = J.yRange;
        M1 = J.marginal1;   M2 = J.marginal2;
    
        if CDF; z = J.cdf; title = "cdf" else z = J.density; title = "pdf";end
    
        m = size(z)[1];
        if m < pn; ppn = m; else ppn = pn; end
    
        nm = round(m/ppn);
        z = z[1:Int(nm):end,1:Int(nm):end]
    
        x = x[1:Int(nm):end]
        y = y[1:Int(nm):end]
    
        xgrid = repeat(x',pn,1)
        ygrid = repeat(y,1,pn)
    
        fig = plt.figure("ContourPlots", figsize=(10, 10))
        grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
        main_ax = fig.add_subplot(get(grid, (slice(1,4),slice(0,3))))
        yMarg  = fig.add_subplot(get(grid, (slice(1,4),3)), xticklabels=[], sharey=main_ax)
        xMarg  = fig.add_subplot(get(grid, (0,slice(0,3))), yticklabels=[], sharex=main_ax)
    
        cp = main_ax.contour(xgrid, ygrid, z, cmap=ColorMap("coolwarm"), levels =15)
        main_ax.clabel(cp, inline=1, fontsize=10)
        #xlabel("X")
        #ylabel("Y")
        if (CDF) xMarg.plot(x,cdf.(M1,x));else xMarg.plot(x,pdf.(M1,x));end
        if (CDF) yMarg.plot(cdf.(M2,y),y);else yMarg.plot(pdf.(M2,y),y);end
    
        PyPlot.title("joint $title and marginals")
        tight_layout()
    
    end



function plotDist(D :: AbstractMarginal; name = "marginal", col = "red", fontsize = 12)

    fig = figure(name,figsize=(10,10))

    PyPlot.plot(D.range, D.cdf, color = col)

    xticks(fontsize = fontsize); yticks(fontsize = fontsize)
    xlabel("Distribution range",fontsize = fontsize); ylabel("CDF",fontsize=fontsize);

end

slice(x,y) = pycall(pybuiltin("slice"), PyObject, x, y)         # For performing python array slicing for

function scatter(a :: Array{Float64,2}; title = "samples",fontsize = 12, fontsizeT = 12, xlab = "U", ylab = "V", nbins = 50)

    x = a[:,1];
    y = a[:,2];

    fig = plt.figure(title, figsize=(10, 10))
    grid = plt.GridSpec(4, 4, hspace=0.2, wspace=0.2)
    main_ax = fig.add_subplot(get(grid, (slice(1,4),slice(0,3))))

    # scatter points on the main axes
    main_ax.plot(x, y, "o", markersize=3, alpha=0.2)

    xticks(fontsize=fontsizeT); yticks(fontsize=fontsizeT)
    xlabel(xlab, fontsize =  fontsize); ylabel(ylab, fontsize = fontsize);

    x_hist  = fig.add_subplot(get(grid, (0,slice(0,3))), xticklabels=[], yticklabels=[])#, sharex=main_ax)

    # histogram on the attached axes
    x_hist.hist(x, nbins, histtype="stepfilled", orientation="vertical", color="gray")
    #x_hist.invert_yaxis()
    xticks(fontsize=fontsizeT); yticks(fontsize=fontsizeT)

    y_hist  = fig.add_subplot(get(grid, (slice(1,4),3)), xticklabels=[], yticklabels=[])#, sharey=main_ax)

    y_hist.hist(y, nbins, histtype="stepfilled", orientation="horizontal", color="gray")
    #y_hist.invert_xaxis()
    xticks(fontsize=fontsizeT); yticks(fontsize=fontsizeT)

end
