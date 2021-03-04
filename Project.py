# -*- coding: utf-8 -*-
"""
Created on Fri Oct 23 03:00:22 2020

@author: Aya Arbiat
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 18:26:32 2020

@author: Aya Arbiat
"""
import math
from sympy import symbols
import numpy as np
from  scipy.interpolate import interp1d
from matplotlib.lines import Line2D
from matplotlib.artist import Artist


class PolygonInteractor(object):

    
    showverts = True
    epsilon = 5  # max pixel distance to count as a vertex hit
    
    def __init__(self, ax, poly, visible=False):
        if poly.figure is None:
            raise RuntimeError('You must first add the polygon to a figure '
                               'or canvas before defining the interactor')
        #CANVAS----------------------------------------------------------------
        self.ax = ax
        canvas = poly.figure.canvas
        self.poly = poly
        self.poly.set_visible(visible)
        
        #POINTS----------------------------------------------------------------
        x, y = zip(*self.poly.xy)
        self.line = Line2D(x, y, ls="",
                           marker='o',markersize=10.5,
                            markerfacecolor='#5c5c5c',animated=True)
        self.ax.add_line(self.line)
    
        self.cid = self.poly.add_callback(self.poly_changed)
        self._ind = None  # the active vert
    
        #CALLS-----------------------------------------------------------------
        canvas.mpl_connect('draw_event', self.draw_callback)
        canvas.mpl_connect('button_press_event', self.button_press_callback)
        canvas.mpl_connect('key_press_event', self.key_press_callback)
        canvas.mpl_connect('button_release_event', self.button_release_callback)
        canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
        self.canvas = canvas
        
        #LINE------------------------------------------------------------------
        x,y = self.Bezier()
        self.line2 = Line2D(x, y,color='#6b6b6b' ,animated=True)
        self.ax.add_line(self.line2)
    
    
    #MATH----------------------------------------------------------------------
    def Bezier(self): 
        
        x, y = self.poly.xy[:].T
        t=symbols('t')
        
        Px =x
        Py =y
        print('len is',len(y))
        n=4
        print(n)
        
        Bi=[]
        Bx=[]
        By=[]
        fx=[]
        fy=[]
        
        for k in [0,1,2,3,4]:
            print(k)
            Bi.append(math.factorial(n) / float(math.factorial(k) * math.factorial(n- k)))
            Bx.append(Px[k]*Bi[k] * (t ** k) * ((1 - t) ** (n- k)))
            By.append(Bi[k] *Py[k]*(t ** k) * ((1 - t) ** (n- k)))
           
        fx=sum(Bx) 
        fy=sum(By) 
        print(fx)
        
        fxs=[]
        fys=[]
        
        step = 0.05
        m= [(x * step) for x in range(0, 21)]
        for j in range(0, len(m)): 
        
            fxs.append(fx.subs(t,m[j]))
            fys.append(fy.subs(t,m[j]))
        
        fxi=np.array(fxs)
        fyi=np.array(fys)  
        
        return fxi, fyi
    
    
    
    
    #GUI-----------------------------------------------------------------------
    def draw_callback(self, event):
        self.background = self.canvas.copy_from_bbox(self.ax.bbox)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.line2)
        # do not need to blit here, this will fire before the screen is
        # updated
    
    def poly_changed(self, poly):
        'this method is called whenever the polygon object is called'
        # only copy the artist props to the line (except visibility)
        vis = self.line.get_visible()
        Artist.update_from(self.line, poly)
        self.line.set_visible(vis)  # don't use the poly visibility state
    
    
    def get_ind_under_point(self, event):
        'get the index of the vertex under point if within epsilon tolerance'
    
        # display coords
        xy = np.asarray(self.poly.xy)
        xyt = self.poly.get_transform().transform(xy)
        xt, yt = xyt[:, 0], xyt[:, 1]
        d = np.hypot(xt - event.x, yt - event.y)
        indseq, = np.nonzero(d == d.min())
        ind = indseq[0]
    
        if d[ind] >= self.epsilon:
            ind = None
    
        return ind
    
    def button_press_callback(self, event):
        'whenever a mouse button is pressed'
        if not self.showverts:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        self._ind = self.get_ind_under_point(event)
    
    def button_release_callback(self, event):
        'whenever a mouse button is released'
        if not self.showverts:
            return
        if event.button != 1:
            return
        self._ind = None
    
    def key_press_callback(self, event):
        'whenever a key is pressed'
        if not event.inaxes:
            return
        if event.key == 't':
            self.showverts = not self.showverts
            self.line.set_visible(self.showverts)
            if not self.showverts:
                self._ind = None
        elif event.key == 'd':
            ind = self.get_ind_under_point(event)
            if ind is not None:
                self.poly.xy = np.delete(self.poly.xy,
                                         ind, axis=0)
                self.line.set_data(zip(*self.poly.xy))
        elif event.key == 'i':
            xys = self.poly.get_transform().transform(self.poly.xy)
            p = event.x, event.y  # display coords
            for i in range(len(xys) - 1):
                s0 = xys[i]
                s1 = xys[i + 1]
                d = dist_point_to_segment(p, s0, s1)
                if d <= self.epsilon:
                    self.poly.xy = np.insert(
                        self.poly.xy, i+1,
                        [event.xdata, event.ydata],
                        axis=0)
                    self.line.set_data(zip(*self.poly.xy))
                    break
        if self.line.stale:
            self.canvas.draw_idle()
    
    def motion_notify_callback(self, event):
        'on mouse movement'
        if not self.showverts:
            return
        if self._ind is None:
            return
        if event.inaxes is None:
            return
        if event.button != 1:
            return
        x, y = event.xdata, event.ydata
    
        self.poly.xy[self._ind] = x, y
        if self._ind == 0:
            self.poly.xy[-1] = x, y
        elif self._ind == len(self.poly.xy) - 1:
            self.poly.xy[0] = x, y
        self.line.set_data(zip(*self.poly.xy))
    
        x,y = self.Bezier()
        self.line2.set_data(x,y)
    
        self.canvas.restore_region(self.background)
        self.ax.draw_artist(self.poly)
        self.ax.draw_artist(self.line)
        self.ax.draw_artist(self.line2)
        self.canvas.blit(self.ax.bbox)



#Callback to class-------------------------------------------------------------
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    
    #INITAL SATATE-------------------------------------------------------------
    xs = (0.5,1.5,3,5,6)
    ys = (0.5,5,2,4,9)

    poly = Polygon(list(zip(xs, ys)),animated=True)
    
    
    fig, ax = plt.subplots()
    ax.add_patch(poly)
    p = PolygonInteractor(ax, poly, visible=False)

    # ax.set_title('Click and drag a point to move it')

    ax.set_xlim((0, 10))
    ax.set_ylim((0, 10))
    plt.grid()
    ax.grid(which='minor', alpha=0.2)
    plt.show()
    
    # major_ticks = np.arange(0, 10, 0.6)
    # minor_ticks = np.arange(0, 10, 0.6)
    
    # ax.set_xticks(major_ticks)
    # ax.set_xticks(minor_ticks, minor=True)
    # ax.set_yticks(major_ticks)
    # ax.set_yticks(minor_ticks, minor=True)


    # ax.grid(which='minor', alpha=0.5)
    # ax.grid(which='major', alpha=0.5)
    
    # plt.show()