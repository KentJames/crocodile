import numpy
import scipy
import scipy.special

from clean import *
from synthesis import *
from simulate import *

from matplotlib import pylab
import matplotlib.pyplot as pyplot

def aaf_ns(a, m, c):
    """

    """
    r=numpy.hypot(*ucs(a))
    return scipy.special.pro_ang1(m,m,c,r)



if 1:
    vlas=numpy.genfromtxt("/home/vlad/software/SKA/crocodile/test/VLA_A_hor_xyz_2ants.txt", delimiter=",")
    pyplot.scatter(vlas[:,0], vlas[:,1])
    pyplot.show()

    vobs=genuv(vlas, numpy.arange(0,numpy.pi,0.1) ,  numpy.pi/4)
    pyplot.scatter(vobs[:,0], vobs[:,1])
    pyplot.show()

    yy=genvis(vobs/5, 0.0, 0.0)
    #yy=yy + genvis(vobs/5, -0.01, -0.01)

# Fill the conjugated part V(-u,-v) = V*(u,v)
    vobs_tmp = numpy.copy(vobs)
    vobs_tmp[:,0] *= -1.0
    vobs_tmp[:,1] *= -1.0
    vobs_new = numpy.concatenate((vobs, vobs_tmp))
    pyplot.scatter(vobs_new[:,0], vobs_new[:,1])
    pyplot.show()


    yy_tmp = numpy.conjugate(yy)
    yy_new = numpy.concatenate((yy, yy_tmp))

    mat_a = numpy.zeros((512,512),'D')
    maxvobs = max(vobs_new[:,0:1])[0] + 1
    nuv = numpy.zeros((512, 512))
    mat_a = grid1(mat_a,(vobs_new/maxvobs),yy_new)
    pyplot.contour(abs(mat_a))
    pyplot.show()

    #mat_a = wgrid(mat_a,(vobs/maxvobs),yy,nuv)


    mat_b = numpy.fft.fftshift(mat_a)
    pyplot.contour(abs(mat_b))
    pyplot.show()


    c = numpy.fft.ifft2(mat_b)
    pyplot.contour(abs(c))
    pyplot.show()

    

if 0:
    majorcycle(0.025, 15000, vobs/5 , yy, 0.1, 5, 100, 250000)
    

if 0: # some other testing code bits
    mg=exmid(numpy.fft.fftshift(numpy.fft.fft2(aaf(a, 0, 3))),5)
    ws=numpy.arange( p[:,2].min(), p[:,2].max(), wstep)
    wr=zip(ws[:-1], ws[1:]) + [ (ws[-1], p[:,2].max() ) ]
    yy=genvis(vobs/5, 0.001, 0.001)
    d,p=doimg(0.025, 15000, vobs/5, yy, lambda *x: wslicimg(*x, wstep=250))
    pylab.matshow(p[740:850,740:850]); pylab.colorbar(); pylab.show()
    x=numpy.zeros_like(d)
    x[1050,1050]=1
    xuv=numpy.fft.fftshift(numpy.fft.fft2(numpy.fft.ifftshift(x)))
