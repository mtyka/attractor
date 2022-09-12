# Attractor.py - Blender add-on that generates various strange attractors
# Copyright (C) 2014 Mike Tyka
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# Corrected by Spirou4D@laposte.net, 2022-12-09.

# 12/09/2022 7bitretro/Googylip added TSUCS1 Attractor

import bpy
import math

from bpy.props import BoolProperty
from bpy.props import FloatProperty
from bpy.props import FloatVectorProperty
from bpy.props import IntProperty

bl_info = {
  "author": "mike.tyka@gmail.com",
  "blender": (2, 90, 0),
  "category": "Add Curve",
  "description": "Creates a strange attractor curve",
  "location": "View3D > Add > Curve > Attractors",
  "name": "Strange Attractors",
  "version": (0, 1),
  "warning": "",
  "wiki_url": ""
}

def get_prop(cls, name, default):
    return bpy.props.FloatProperty(attr=cls + "_" + name,
                                 name=name,
                                 description="",
                                 default=default)

def get_int_prop(cls, name, default):
    return bpy.props.IntProperty(attr=cls + "_" + name,
                                 name=name,
                                 description="",
                                 default=default)
def get_npoints(default=10000):
    return IntProperty(attr="npoints",
                     name="Vertices",
                     description="",
                     min=1, 
                     soft_min=1, 
                     default=default)

class Object_OT_Attractor(bpy.types.Operator):
    bl_idname= "object.attractor"
    bl_label= "Label Attractor"
    view_align : BoolProperty(name="Align to View",default=False)
    location : FloatVectorProperty(name="Location",subtype="TRANSLATION")
    rotation : FloatVectorProperty(name="Rotation",subtype="EULER")

    def draw(self, context):
        layout = self.layout
        row = layout.row()
        row.prop(self, "npoints", text="Vertices")

        col = layout.column(align=True)
        for param in self.params:
            col.prop(self, param, text=param)

        col2 = layout.column(align=True)
        col2.prop(self, "x", text="x")
        col2.prop(self, "y", text="y")
        col2.prop(self, "z", text="z")
        col2.prop(self, "dt", text="dt")


    def execute(self, context):
        curvedata = bpy.data.curves.new(name="Curve", type="CURVE")
        curvedata.dimensions = "3D"

        objectdata = bpy.data.objects.new("ObjCurve", curvedata)
        objectdata.location = (0,0,0)
        bpy.context.scene.collection.objects.link(objectdata)

        polyline = curvedata.splines.new("POLY")
        
        polyline.points.add(self.npoints - 1)  

        x = self.x
        y = self.y
        z = self.z

        n = 0
        while n < self.npoints:
            dx, dy, dz = self.iterate(x, y, z)
            x = x + self.dt * dx
            y = y + self.dt * dy
            z = z + self.dt * dz
            polyline.points[n].co = (x,y,z,math.sqrt(dx*dx+dy*dy+dz*dz))
            n += 1

        return {"FINISHED"}

class CoulletAttractor(Object_OT_Attractor):
    bl_idname = "curve.coullet_attractor_add"
    bl_label = "Coullet Attractor"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c", "d"]
    a : get_prop("coullet", "a",  0.80)
    b : get_prop("coullet", "b", -1.10)
    c : get_prop("coullet", "c", -0.45)
    d : get_prop("coullet", "d", -1.00)

    x : get_prop("coullet", "x",   0.10)
    y : get_prop("coullet", "y",   0.00)
    z : get_prop("coullet", "z",   0.00)
    dt : get_prop("coullet", "dt", 0.01)

    def iterate(self, x, y, z):
        dx = y
        dy = z
        dz = self.a*x + self.b*y + self.c*z + self.d*x**3
        return (dx, dy, dz)

class LorenzAttractor(Object_OT_Attractor):
    bl_idname = "curve.lorenz_attractor_add"
    bl_label = "Lorenz Attractor"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c"]
    a : get_prop("lorenz", "a",  10.00)
    b : get_prop("lorenz", "b",  28.00)
    c : get_prop("lorenz", "c",  (8/3))

    x : get_prop("lorenz", "x",   0.10)
    y : get_prop("lorenz", "y",   0.00)
    z : get_prop("lorenz", "z",   0.00)
    dt : get_prop("lorenz", "dt", 0.01)

    def iterate(self, x, y, z):
        dx = self.a * (y - x)
        dy = (x * (self.b - z) - y)
        dz = (x * y - self.c * z)
        return (dx, dy, dz)

class RoesslerAttractor(Object_OT_Attractor):
    ## Herrmann, Dietmar: Algorithmen für Fraktale und Chaostheorie
    ## Addison Wesley, 1994 ISBN 3-89319-633-1
    bl_idname = "curve.roessler_attractor_add"
    bl_label = "RoesslerAttractor"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c"]
    a : get_prop("roessler", "a", 0.2)
    b : get_prop("roessler", "b", 0.2)
    c : get_prop("roessler", "c", 5.7)

    x : get_prop("roessler", "x",   1.00)
    y : get_prop("roessler", "y",   1.00)
    z : get_prop("roessler", "z",   1.00)
    dt : get_prop("roessler", "dt", 0.01)

    def iterate(self, x, y, z):
        dx=-(y+z);
        dy=x+self.a*y;
        dz=self.b+z*(x-self.c);
        return (dx, dy, dz)

class AizawaAttractor(Object_OT_Attractor):
    bl_idname = "curve.aizawa_attractor_add"
    bl_label = "Aizawa Attractor"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c", "d", "e", "f"]
    a : get_prop("aizawa", "a",  0.95)
    b : get_prop("aizawa", "b",  0.7 )
    c : get_prop("aizawa", "c",  0.6 )
    d : get_prop("aizawa", "d",  3.5 )
    e : get_prop("aizawa", "e",  0.25)
    f : get_prop("aizawa", "f",  0.1 )
    x : get_prop("aizawa", "x",  0.1 )
    y : get_prop("aizawa", "y",  0.00)
    z : get_prop("aizawa", "z",  0.00)
    dt : get_prop("aizawa", "dt", 0.01)

    def iterate(self, x, y, z):
        dx = (z-self.b)*x-self.d*y
        dy = self.d*x+(z-self.b)*y
        dz = self.c+self.a*z-(pow(z,3)/3)-(pow(x,2)+pow(y,2))*(1+self.e*z)+self.f*z*pow(x,3)
        return (dx, dy, dz)

class ActAttractor(Object_OT_Attractor):
    bl_idname = "curve.act_attractor_add"
    bl_label = "act"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "d", "m"]
    a : get_prop("act", "a", 1.8)
    b : get_prop("act", "b", -0.07)
    d : get_prop("act", "d", 1.5)
    m : get_prop("act", "m", 0.02)
    x : get_prop("act", "x", 0.5)
    y : get_prop("act", "y", 0.0)
    z : get_prop("act", "z", 0.0)
    dt : get_prop("act", "dt", 0.02)

    def iterate(self, x, y, z):
        xn = self.a*(x-y)
        yn = -4*self.a*y+x*z+self.m*pow(x,3)
        zn = -self.d*self.a*z+x*y+self.b*z*z
        return xn, yn, zn

class ThreeCellsCNNAttractor(Object_OT_Attractor):
    bl_idname = "curve.threecellscnn_attractor_add"
    bl_label = "ThreeCellsCNN Attractor"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["p1", "p2", "r", "s"]
    p1 : get_prop("ThreeCellsCNN", "p1", 1.24)
    p2 : get_prop("ThreeCellsCNN", "p2", 1.1)
    r : get_prop("ThreeCellsCNN", "r", 4.4)
    s : get_prop("ThreeCellsCNN", "s", 3.21)

    x : get_prop("ThreeCellsCNN", "x", 0.1)
    y : get_prop("ThreeCellsCNN", "y", 0.1)
    z : get_prop("ThreeCellsCNN", "z", 0.1)
    dt : get_prop("ThreeCellsCNN", "dt", 0.01)

    def iterate(self, x, y, z):
        h1 = 0.5*(abs(x+1)-abs(x-1))
        h2 = 0.5*(abs(y+1)-abs(y-1))
        h3 = 0.5*(abs(z+1)-abs(z-1))
        xn = -x+self.p1*h1-self.s*h2-self.s*h3
        yn = -y-self.s*h1+self.p2*h2-self.r*h3
        zn = -z-self.s*h1+self.r*h2+h3
        return xn, yn, zn

class ArneodoAttractor(Object_OT_Attractor):
    bl_idname = "curve.arneodo_attractor_add"
    bl_label = "Arneodo"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c"]
    a : get_prop("Arneodo", "a", -5.5)
    b : get_prop("Arneodo", "b", 3.5)
    c : get_prop("Arneodo", "c", -1)
    x : get_prop("Arneodo", "x", 0.1)
    y : get_prop("Arneodo", "y", 0.0)
    z : get_prop("Arneodo", "z", 0.0)
    dt : get_prop("Arneodo", "dt", 0.009)

    def iterate(self, x, y, z):
        xn = y
        yn = z
        zn = -self.a*x-self.b*y-z+self.c*pow(x,3)
        return xn, yn, zn

class AnishenkoAstakhovAttractor(Object_OT_Attractor):
    bl_idname = "curve.anishenko_astakhov_attractor_add"
    bl_label = "AnishenkoAstakhov"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["m", "n"]
    m : get_prop("AnishenkoAstakov", "m", 1.2)
    n : get_prop("AnishenkoAstakov", "n", 0.5)
    x : get_prop("AnishenkoAstakov", "x", 1)
    y : get_prop("AnishenkoAstakov", "y", 0.0)
    z : get_prop("AnishenkoAstakov", "z", 0.0)
    dt : get_prop("AnishenkoAstakov", "dt", 0.01)

    def iterate(self, x, y, z):
        i = 1.0 if x > 0 else 0.0
        xn = self.m*x + y - x*z
        yn = -x
        zn = -self.n*z + self.n*i*x*x
        return xn, yn, zn

class BoualiAttractor(Object_OT_Attractor):
    bl_idname = "curve.bouali_attractor_add"
    bl_label = "Bouali"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "g"]
    a : get_prop("Bouali", "a", 0.2)
    b : get_prop("Bouali", "b", 0.05)
    g : get_prop("Bouali", "g", 1.0)

    x : get_prop("Bouali", "x", 1)
    y : get_prop("Bouali", "y", 0.0)
    z : get_prop("Bouali", "z", 0.0)
    dt : get_prop("Bouali", "dt", 0.005)

    def iterate(self, x, y, z):
        xn = 0.02*y + 0.4*x*(0.2 - y*y)
        yn = -x + self.a*z
        zn = 10.0*x - 0.1*y
        return xn, yn, zn

class BurkeShawAttractor(Object_OT_Attractor):
    bl_idname = "curve.burkeshaw_attractor_add"
    bl_label = "BurkeShaw"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["s", "v"]
    s : get_prop("BurkeShaw", "s", 10)
    v : get_prop("BurkeShaw", "v", 4.272)
    x : get_prop("BurkeShaw", "x", 1)
    y : get_prop("BurkeShaw", "y", 0.0)
    z : get_prop("BurkeShaw", "z", 0.0)
    dt : get_prop("BurkeShaw", "dt", 0.005)

    def iterate(self, x, y, z):
        xn = -self.s*(x+y)
        yn = -y-self.s*x*z
        zn = self.s*x*y+self.v
        return xn, yn, zn

class ChenAttractor(Object_OT_Attractor):
    ## Professor Tetsushi Ueta
    ## Center for Advanced Information Technology, Tokushima University
    bl_idname = "curve.chen_attractor_add"
    bl_label = "Chen"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c"]
    a : get_prop("Chen", "a", 35)
    b : get_prop("Chen", "b", 8/3.0)
    c : get_prop("Chen", "c", 20)
    x : get_prop("Chen", "x", -3)
    y : get_prop("Chen", "y", 2)
    z : get_prop("Chen", "z", 20)
    dt : get_prop("Chen", "dt", 0.002)

    def iterate(self, x, y, z):
        xn = self.a*(y-x)
        yn = (self.c-self.a)*x-(x*z)+(self.c*y)
        zn = (x*y)-(self.b*z)
        return xn, yn, zn

class LotkaVolterraAttractor(Object_OT_Attractor):
    ## John S. Costello
    ## Synchronisation of Chaos in a Generalized Lotka-Volterra Atractor
    ## The Nonlinear Journal Vol. 1, 1999 pp 11-17
    bl_idname = "curve.lotka_volterra_attractor_add"
    bl_label = "LotkaVolterra"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints(30000)
    params = ["a", "b", "c"]
    a : get_prop("LotkaVolterra", "a", 2.9851)
    b : get_prop("LotkaVolterra", "b", 3.0)
    c : get_prop("LotkaVolterra", "c", 2)
    x : get_prop("LotkaVolterra", "x", 1.0)
    y : get_prop("LotkaVolterra", "y", 1.0)
    z : get_prop("LotkaVolterra", "z", 1.0)
    dt : get_prop("LotkaVolterra", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = x - x*y + self.c*x*x - self.a*z*x*x
        yn = -y + x*y
        zn = -self.b*z + self.a*z*x*x
        return xn, yn, zn

class MooreSpiegelAttractor(Object_OT_Attractor):
    ## L. A. Smith, C. Ziehmann und K. Fraedrich
    ## Uncertainty Dynamics and Predictability in Chaotic Systems
    ## (1999) Quart. J. Royal Meteorological Soc. 125, 2855-2886
    ## http://www.lse.ac.uk/collections/cats/abstracts/UncertaintyDynamics.htm
    bl_idname = "curve.moorespiegel_attractor_add"
    bl_label = "MooreSpiegel"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b"]
    a : get_prop("MoreSpiegel", "a", 100.0)
    b : get_prop("MoreSpiegel", "b", 26.0)

    x : get_prop("MoreSpiegel", "x", 1.0)
    y : get_prop("MoreSpiegel", "y", 0.0)
    z : get_prop("MoreSpiegel", "z", 0.0)
    dt : get_prop("MoreSpiegel", "dt", 0.002)

    def iterate(self, x, y, z):
        xn=y
        yn=z
        zn=(-z-(self.b-self.a+self.a*x*x)*y-self.b*x)
        return xn, yn, zn

class RikitakeAttractor(Object_OT_Attractor):
    ## Tyler McMillen
    ## The Shape and Dynamics of the Rikitake Atractor
    ## The Nonlinear Journal Vol. 1, 1999 pp 1-10
    bl_idname = "curve.rikitake_attractor_add"
    bl_label = "Rikitake"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b"]
    a : get_prop("Rikitake", "a", 5.0)
    b : get_prop("Rikitake", "b", 2.0)
    x : get_prop("Rikitake", "x", 1.0)
    y : get_prop("Rikitake", "y", 1.0)
    z : get_prop("Rikitake", "z", 1.0)
    dt : get_prop("Rikitake", "dt", 0.01)

    def iterate(self, x, y, z):
        xn=-self.b*x+z*y
        yn=-self.b*y+(z-self.a)*x
        zn=1.0-x*y
        return xn, yn, zn

class RabinovichFabrikantAttractor(Object_OT_Attractor):
    ## Rabinovich Fabrikant Attraktor: Robert Doerner
    ## http://www.robert-doerner.de/Glossar/glossar.html#rabinovich
    bl_idname = "curve.rabinovichfabrikant_attractor_add"
    bl_label = "RabinovichFabrikant"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints(25000)
    params = ["a", "b"]
    a : get_prop("RabinovichFabrikant", "a", 1.1)
    b : get_prop("RabinovichFabrikant", "b", 0.87)
    x : get_prop("RabinovichFabrikant", "x", -1.05)
    y : get_prop("RabinovichFabrikant", "y", 0.9)
    z : get_prop("RabinovichFabrikant", "z", 1.01)
    dt : get_prop("RabinovichFabrikant", "dt", 0.001)

    def iterate(self, x, y, z):
        xn=y*(z-1.0+x*x)+self.b*x
        yn=x*(3.0*z+1.0-x*x)+self.b*y
        zn=-2.0*z*(self.a+x*y)
        return xn, yn, zn

class ThreeLayerAttractor(Object_OT_Attractor):
    ## Tianshou Zhou, Guanrong Chen A Simple Smooth Chaotic System with a 3-Layer Attractor
    bl_idname = "curve.threelayer_attractor_add"
    bl_label = "ThreeLayer"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints(15000)
    params = ["a1","a2","a3","c1","c2","c3","d"]
    a1 : get_prop("ThreeLayer", "a1",-4.1)
    a2 : get_prop("ThreeLayer", "a2", 1.2)
    a3 : get_prop("ThreeLayer", "a3",13.45)
    c1 : get_prop("ThreeLayer", "c1", 2.76)
    c2 : get_prop("ThreeLayer", "c2", 0.6)
    c3 : get_prop("ThreeLayer", "c3", 13.13)
    d : get_prop("ThreeLayer", "d", 1.8)

    x : get_prop("ThreeLayer", "x", 0.0)
    y : get_prop("ThreeLayer", "y", -15.0)
    z : get_prop("ThreeLayer", "z", 0.0)
    dt : get_prop("ThreeLayer", "dt", 0.04 )

    def iterate(self, x, y, z):
        b=((self.d*self.a2*self.a2*self.c3*self.c3)/(32*self.a3*self.a3*self.c2*self.c2))*math.sqrt((-self.a3*self.c2)/(self.a1*self.c1))
        c=((self.a2*self.a2*self.c3*self.c3)/(4*self.a3*self.c2) - (self.a3*self.c1 + self.a1*self.c2)*b/self.d)/self.a2
        xn=self.a1*x - self.a2*y + self.a3*z
        yn=-self.d*x*z + b
        zn=self.c1*x*y + self.c2*y*z + self.c3*z + c
        return xn, yn, zn

class ChuaAttractor(Object_OT_Attractor):
    ## http://haides.caltech.edu/~mcc/chaos_new/Chua.html
    ## http://www.science.gmu.edu/~rcastro/chualink.html
    bl_idname = "curve.chua_attractor_add"
    bl_label = "Chua"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a","b","c","d","e"]
    a : get_prop("Chua", "a", 15.6)
    b : get_prop("Chua", "b", 1.0)
    c : get_prop("Chua", "c", 25.58)
    d : get_prop("Chua", "d", -1)
    e : get_prop("Chua", "e", 0)

    x : get_prop("Chua", "x", 1.0)
    y : get_prop("Chua", "y", 1.0)
    z : get_prop("Chua", "z", 1.0)
    dt : get_prop("Chua", "dt", 0.01)

    def iterate(self, x, y, z):
        h = self.e*x+(self.d+self.e)*(abs(x+1)-abs(x-1));
        xn = self.a * (y - x - h);
        yn = self.b * (x - y + z);
        zn = -self.c * y;
        return xn, yn, zn

class ChuaMultiIAttractor(Object_OT_Attractor):
    ## W. Tang, G. Zhong, G. Chen und K. Man Generation of N-Scroll Attractors via Sine Function
    bl_idname = "curve.chuamultii_attractor_add"
    bl_label = "ChuaMultiI"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints(15000)
    params = ["alpha","beta","a","b","c"]
    alpha : get_prop("ChuaMultiI", "alpha", 10.814)
    beta : get_prop("ChuaMultiI", "beta", 14.0)
    a : get_prop("ChuaMultiI", "a", 1.3)
    b : get_prop("ChuaMultiI", "b", 0.11)
    c : get_int_prop("ChuaMultiI", "c", 2)

    x : get_prop("ChuaMultiI", "x", 1.0)
    y : get_prop("ChuaMultiI", "y", 0.0)
    z : get_prop("ChuaMultiI", "z", 0.0)
    dt : get_prop("ChuaMultiI", "dt", 0.01)

    def iterate(self, x, y, z):
        d=0
        if self.c % 2 == 0: d=math.pi
        h = 0
        if x >= 2*self.a*self.c:
             h=((self.b*math.pi)/(2*self.a))*(x - 2*self.a*self.c)
        if x <= -2*self.a*self.c:
             h=((self.b*math.pi)/(2*self.a))*(x + 2*self.a*self.c)
        if x < 2*self.a*self.c  and  x > -2*self.a*self.c:
             h=-self.b*math.sin((math.pi*x)/(2*self.a)+d)

        xn = self.alpha * (y - h)
        yn = x - y + z
        zn = -self.beta * y
        return xn, yn, zn

class ChuaMultiIIAttractor(Object_OT_Attractor):
    ## J.A.K. Suykens und J. Vandewalle
    ## The K.U. Leuven Competition Data: a Challenge for Advanced Neural Network Techniques
    ## http://www.dice.ucl.ac.be/Proceedings/esann/esannpdf/es2000-261.pdf
    bl_idname = "curve.chuamultiii_attractor_add"
    bl_label = "ChuaMultiII"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints(15000)
    params = ["alpha","beta","a","b"]
    alpha : get_prop("ChuaMultiII", "alpha", 10.814)
    beta : get_prop("ChuaMultiII", "beta", 14.0)
    a : get_prop("ChuaMultiII", "a", 1.3)
    b : get_prop("ChuaMultiII", "b", 0.11)

    x : get_prop("ChuaMultiII", "x", 1.0)
    y : get_prop("ChuaMultiII", "y", 0.0)
    z : get_prop("ChuaMultiII", "z", 0.0)
    dt : get_prop("ChuaMultiII", "dt", 0.01)

    m = [0.9/7.0, -3.0/7.0, 3.5/7.0, -2.7/7.0, 4.0/7.0, -2.4/7.0]
    c = [0.0, 1.0, 2.15, 3.6, 6.2, 9.0]

    def iterate(self, x, y, z):
        sum=0
        for k in range(1,6):
             sum+=((self.m[k-1] - self.m[k])*(abs(x + self.c[k])-abs(x - self.c[k])))
        h = self.m[5]*x+0.5*sum
        xn = self.alpha * (y - h)
        yn = x - y + z
        zn = -self.beta * y
        return xn, yn, zn

class ChenLeeAttractor(Object_OT_Attractor):
    bl_idname = "curve.chenlee_attractor_add"
    bl_label = "ChenLee"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c"]
    a : get_prop("ChenLee", "a", 5)
    b : get_prop("ChenLee", "b", -10)
    c : get_prop("ChenLee", "c", -0.38)
    x : get_prop("ChenLee", "x", 1)
    y : get_prop("ChenLee", "y", 1)
    z : get_prop("ChenLee", "z", 1)
    dt : get_prop("ChenLee", "dt", 0.002)

    def iterate(self, x, y, z):
        xn = self.a*x - y*z
        yn = self.b*y + x*z
        zn = self.c*z + x*y/3.0
        return xn, yn, zn

class FinanceAttractor(Object_OT_Attractor):
    bl_idname = "curve.finance_attractor_add"
    bl_label = "Finance"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c"]
    a : get_prop("Finance", "a", 0.001)
    b : get_prop("Finance", "b", 0.2)
    c : get_prop("Finance", "c", 1.1)
    x : get_prop("Finance", "x", 0.1)
    y : get_prop("Finance", "y", 0.1)
    z : get_prop("Finance", "z", 0.1)
    dt : get_prop("Finance", "dt", 0.002)

    def iterate(self, x, y, z):
        xn = (1/self.b - self.a)*x + x*y + z
        yn = -self.b*y - x*x
        zn = -x-self.c*z
        return xn, yn, zn


class ChenCelikovskyAttractor(Object_OT_Attractor):
    bl_idname = "curve.chencelikovsky_attractor_add"
    bl_label = "ChenCelikovsky"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c"]
    a : get_prop("ChenCelikovsky", "a", 36)
    b : get_prop("ChenCelikovsky", "b", 3)
    c : get_prop("ChenCelikovsky", "c", 20)
    x : get_prop("ChenCelikovsky", "x", 1)
    y : get_prop("ChenCelikovsky", "y", 1)
    z : get_prop("ChenCelikovsky", "z", 1)
    dt : get_prop("ChenCelikovsky", "dt", 0.002)
    cnt = 0
    t = 0

    def iterate(self, x, y, z):
        xn = self.a*(y-x)
        yn = -(x*z)+(self.c*y)
        zn = (x*y)-(self.b*z)
        self.cnt += 1
        if self.cnt >= 500: self.t = 1
        return xn, yn, zn

class DadrasAttractor(Object_OT_Attractor):
    bl_idname = "curve.dadras_attractor_add"
    bl_label = "Dadras"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["p", "o", "r", "c", "e"]
    p : get_prop("Dadras", "p", 3)
    o : get_prop("Dadras", "o", 2.7)
    r : get_prop("Dadras", "r", 1.7)
    c : get_prop("Dadras", "c", 2)
    e : get_prop("Dadras", "e", 9)
    x : get_prop("Dadras", "x", 0.1)
    y : get_prop("Dadras", "y", 0.03)
    z : get_prop("Dadras", "z", 0.0)
    dt : get_prop("Dadras", "dt", 0.008)

    def iterate(self, x, y, z):
        xn = y - self.p*x + self.o*y*z
        yn = self.r*y - x*z + z
        zn = self.c*x*y - self.e*z
        return xn, yn, zn

class DequanLiAttractor(Object_OT_Attractor):
    bl_idname = "curve.dequanli_attractor_add"
    bl_label = "DequanLi"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints(20000)
    params = ["a", "b", "d", "e", "p", "g"]
    a : get_prop("DequanLi", "a", 40)
    b : get_prop("DequanLi", "b", 11.0/6.0)
    d : get_prop("DequanLi", "d", 0.16)
    e : get_prop("DequanLi", "e", 0.65)
    p : get_prop("DequanLi", "p", 55.0)
    g : get_prop("DequanLi", "g", 20.0)

    x : get_prop("DequanLi", "x", 0.01)
    y : get_prop("DequanLi", "y", 0.0)
    z : get_prop("DequanLi", "z", 0.0)
    dt : get_prop("DequanLi", "dt", 0.0003)

    def iterate(self, x, y, z):
        xn = self.a*(y-x) + self.d*x*z
        yn = self.p*x + self.g*y - x*z
        zn = self.b*z + x*y - self.e*x*x
        return xn, yn, zn

class YuWangAttractor(Object_OT_Attractor):
    bl_idname = "curve.yuwang_attractor_add"
    bl_label = "YuWang"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c", "d"]
    a : get_prop("YuWang", "a", 10.0)
    b : get_prop("YuWang", "b", 40.0)
    c : get_prop("YuWang", "c", 2)
    d : get_prop("YuWang", "d", 2.5)
    x : get_prop("YuWang", "x", 0.1)
    y : get_prop("YuWang", "y", 0.0)
    z : get_prop("YuWang", "z", 15)
    dt : get_prop("YuWang", "dt", 0.002)

    def iterate(self, x, y, z):
        xn = self.a*(y-x)
        yn = self.b * x - self.c*x*z
        zn = math.exp(x*y) - self.d*z
        return xn, yn, zn

class HadleyAttractor(Object_OT_Attractor):
    bl_idname = "curve.hadley_attractor_add"
    bl_label = "Hadley"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "f", "g"]
    a : get_prop("Hadley", "a", 0.2)
    b : get_prop("Hadley", "b", 4)
    f : get_prop("Hadley", "f", 8)
    g : get_prop("Hadley", "g", 1)
    x : get_prop("Hadley", "x", 0.0)
    y : get_prop("Hadley", "y", 0.0)
    z : get_prop("Hadley", "z", 1.3)
    dt : get_prop("Hadley", "dt", 0.008)

    def iterate(self, x, y, z):
        xn = -y*y-z*z-self.a*x+self.a*self.f
        yn = x*y-self.b*x*z-y+self.g
        zn = self.b*x*y+x*z-z
        return xn, yn, zn

class HalvorsenAttractor(Object_OT_Attractor):
    bl_idname = "curve.halvorsen_attractor_add"
    bl_label = "Halvorsen"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    param = ["a"]
    a : get_prop("Halvorsen", "a", 1.4)
    x : get_prop("Halvorsen", "x", -5)
    y : get_prop("Halvorsen", "y", 0.0)
    z : get_prop("Halvorsen", "z", 0.0)
    dt : get_prop("Halvorsen", "dt", 0.004)

    def iterate(self, x, y, z):
        xn = -self.a*x-4*y-4*z-y*y
        yn = -self.a*y-4*z-4*x-z*z
        zn = -self.a*z-4*x-4*y-x*x
        return xn, yn, zn

class LinzSprottAttractor(Object_OT_Attractor):
    bl_idname = "curve.linzsprott_attractor_add"
    bl_label = "LinzSprott"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    param = []
    x : get_prop("LinzSprott", "x", 0.1)
    y : get_prop("LinzSprott", "y", 0.1)
    z : get_prop("LinzSprott", "z", 0.1)
    dt : get_prop("LinzSprott", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = y
        yn = -x+y*z
        zn = 1-pow(y,2)
        return xn, yn, zn

class LorenzMod1Attractor(Object_OT_Attractor):
  bl_idname = "curve.lorenzmod1_attractor_add"
  bl_label = "LorenzMod1"
  bl_options = {"REGISTER", "UNDO"}

  npoints : get_npoints()
  params = ["a", "b", "c", "d"]
  a : get_prop("LorenzMod1", "a", 0.1)
  b : get_prop("LorenzMod1", "b", 4)
  c : get_prop("LorenzMod1", "c", 14)
  d : get_prop("LorenzMod1", "d", 0.08)
  x : get_prop("LorenzMod1", "x", 0.0)
  y : get_prop("LorenzMod1", "y", 1)
  z : get_prop("LorenzMod1", "z", 0.0)
  dt : get_prop("LorenzMod1", "dt", 0.005)

  def iterate(self, x, y, z):
    xn = -self.a*x+(y*y)-(z*z)+self.a*self.c
    yn = x*(y-self.b*z)+self.d
    zn = -z+x*(self.b*y+z)
    return xn, yn, zn

class LorenzMod2Attractor(Object_OT_Attractor):
    bl_idname = "curve.lorenzmod2_attractor_add"
    bl_label = "LorenzMod2"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "c", "d"]
    a : get_prop("LorenzMod2", "a", 0.9)
    b : get_prop("LorenzMod2", "b", 5)
    c : get_prop("LorenzMod2", "c", 9.9)
    d : get_prop("LorenzMod2", "d", 1)
    x : get_prop("LorenzMod2", "x", 5)
    y : get_prop("LorenzMod2", "y", 5)
    z : get_prop("LorenzMod2", "z", 5)
    dt : get_prop("LorenzMod2", "dt", 0.003)

    def iterate(self, x, y, z):
        xn = -self.a*x+(y*y)-(z*z)+self.a*self.c
        yn = x*(y-self.b*z)+self.d
        zn = -z+x*(self.b*y+z)
        return xn, yn, zn

class NewtonLeipnikAttractor(Object_OT_Attractor):
    bl_idname = "curve.newtonleipnik_attractor_add"
    bl_label = "NewtonLeipnik"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b"]
    a : get_prop("NewtonLeipnik", "a", 0.4)
    b : get_prop("NewtonLeipnik", "b", 0.175)
    x : get_prop("NewtonLeipnik", "x", 0.349)
    y : get_prop("NewtonLeipnik", "y", 0.0)
    z : get_prop("NewtonLeipnik", "z", -0.16)
    dt : get_prop("NewtonLeipnik", "dt", 0.03)

    def iterate(self, x, y, z):
        xn = -self.a*x+y+10*y*z
        yn = -x-0.4*y+5*x*z
        zn = self.b*z-5*x*y
        return xn, yn, zn

class NoseHooverAttractor(Object_OT_Attractor):
    bl_idname = "curve.nosehoover_attractor_add"
    bl_label = "NoseHoover"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("NoseHoover", "a", 1.5)
    x : get_prop("NoseHoover", "x", 1)
    y : get_prop("NoseHoover", "y", 0.0)
    z : get_prop("NoseHoover", "z", 0.0)
    dt : get_prop("NoseHoover", "dt", 0.009)

    def iterate(self, x, y, z):
        xn = y
        yn = -x+y*z
        zn = self.a-(y*y)
        return xn, yn, zn

class RayleighBenardAttractor(Object_OT_Attractor):
    bl_idname = "curve.rayleighbenard_attractor_add"
    bl_label = "RayleighBenard"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "r", "b"]
    a : get_prop("RayleighBenard", "a", 9)
    r : get_prop("RayleighBenard", "r", 12)
    b : get_prop("RayleighBenard", "b", 5)
    x : get_prop("RayleighBenard", "x", 0.1)
    y : get_prop("RayleighBenard", "y", 0.0)
    z : get_prop("RayleighBenard", "z", 0.0)
    dt : get_prop("RayleighBenard", "dt", 0.05)

    def iterate(self, x, y, z):
        xn = -self.a*x+self.a*y
        yn = self.r*x-y-x*z
        zn = x*y-self.b*z
        return xn, yn, zn

class RucklidgeAttractor(Object_OT_Attractor):
    bl_idname = "curve.rucklidge_attractor_add"
    bl_label = "Rucklidge"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["k", "a"]
    k : get_prop("Rucklidge", "k", 2)
    a : get_prop("Rucklidge", "a", 6.7)
    x : get_prop("Rucklidge", "x", -5)
    y : get_prop("Rucklidge", "y", 0.0)
    z : get_prop("Rucklidge", "z", 0.0)
    dt : get_prop("Rucklidge", "dt", 0.008)

    def iterate(self, x, y, z):
        xn = -self.k*x+self.a*y-y*z
        yn = x
        zn = -z+y*y
        return xn, yn, zn

class WangSunAttractor(Object_OT_Attractor):
    bl_idname = "curve.wangsun_attractor_add"
    bl_label = "WangSun"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints(30000)
    params = ["a", "b", "c", "d", "e", "f"]
    a : get_prop("WangSun", "a", 0.2)
    b : get_prop("WangSun", "b", -0.01)
    c : get_prop("WangSun", "c", 1.0)
    d : get_prop("WangSun", "d", -0.4)
    e : get_prop("WangSun", "e", -1.0)
    f : get_prop("WangSun", "f", -1.0)
    x : get_prop("WangSun", "x", 0.5)
    y : get_prop("WangSun", "y", 0.1)
    z : get_prop("WangSun", "z", 0.1)
    dt : get_prop("WangSun", "dt", 0.02)

    def iterate(self, x, y, z):
        xn = self.a*x + self.c*y*z
        yn = self.b*x + self.d*y - x*z
        zn = self.e*z + self.f*x*y
        return xn, yn, zn

class NewJerkAttractor(Object_OT_Attractor):
    bl_idname = "curve.newjerk_attractor_add"
    bl_label = "NewJerk"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["k", "a"]
    a : get_prop("NewJerk", "a", 0.00000001)
    b : get_prop("NewJerk", "b", 0.026)
    x : get_prop("NewJerk", "x", 0.1)
    y : get_prop("NewJerk", "y", 0.1)
    z : get_prop("NewJerk", "z", 0.1)
    dt : get_prop("NewJerk", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = y
        yn = z
        zn = -z-x-self.a*(math.exp(y/self.b)-1.0)
        return xn, yn, zn

class SakaryaAttractor(Object_OT_Attractor):
    bl_idname = "curve.sakarya_attractor_add"
    bl_label = "Sakarya"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b"]
    a : get_prop("Sakarya", "a", 0.4)
    b : get_prop("Sakarya", "b", 0.3)
    x : get_prop("Sakarya", "x", 1)
    y : get_prop("Sakarya", "y", -1)
    z : get_prop("Sakarya", "z", 1)
    dt : get_prop("Sakarya", "dt", 0.003)

    def iterate(self, x, y, z):
        xn = -x+y+(y*z)
        yn = -x-y+self.a*(x*z)
        zn = z-self.b*(x*y)
        return xn, yn, zn

class Duffing(Object_OT_Attractor):
    bl_idname = "curve.duffing_attractor_add"
    bl_label = "Duffing"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "w"]
    a : get_prop("Duffing", "a", 0.25)
    b : get_prop("Duffing", "b", 0.3)
    w : get_prop("Duffing", "w", 1.0)
    x : get_prop("Duffing", "x", 0.1)
    y : get_prop("Duffing", "y", 0.0)
    z : get_prop("Duffing", "z", 0.0)
    dt : get_prop("Duffing", "dt", 0.01)
    t : 0.0

    def iterate(self, x, y, z):
        xn = y
        yn = x-x**3-self.a*y+self.b*math.cos(self.t*self.w)
        zn = 0
        t += dt
        return xn, yn, zn

class ShimizuMoriokaAttractor(Object_OT_Attractor):
    bl_idname = "curve.shimizumorioka_attractor_add"
    bl_label = "ShimizuMorioka"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b"]
    a : get_prop("ShimizuMorioka", "a", 0.75)
    b : get_prop("ShimizuMorioka", "b", 0.45)
    x : get_prop("ShimizuMorioka", "x", 0.1)
    y : get_prop("ShimizuMorioka", "y", 0.0)
    z : get_prop("ShimizuMorioka", "z", 0.0)
    dt : get_prop("ShimizuMorioka", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = y
        yn = (1-z)*x-self.a*y
        zn = (x*x)-self.b*z
        return xn, yn, zn

class SprottDAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottd_attractor_add"
    bl_label = "SprottD"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottD", "a", 3)
    x : get_prop("SprottD", "x", 0.1)
    y : get_prop("SprottD", "y", 0.0)
    z : get_prop("SprottD", "z", 0.0)
    dt : get_prop("SprottD", "dt", 0.008)

    def iterate(self, x, y, z):
        xn = -y
        yn = x+z
        zn = x*z+self.a*pow(y,2)
        return xn, yn, zn

class SprottEAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprotte_attractor_add"
    bl_label = "SprottE"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottE", "a", 4)
    x : get_prop("SprottE", "x", 1)
    y : get_prop("SprottE", "y", 0.0)
    z : get_prop("SprottE", "z", 0.0)
    dt : get_prop("SprottE", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = y*z
        yn = pow(x,2)-y
        zn = 1-self.a*x
        return xn, yn, zn

class SprottLinzBAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzb_attractor_add"
    bl_label = "SprottLinzB"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    param = []
    x : get_prop("SprottLinzB", "x", 1)
    y : get_prop("SprottLinzB", "y", 0.0)
    z : get_prop("SprottLinzB", "z", 1)
    dt : get_prop("SprottLinzB", "dt", 0.02)

    def iterate(self, x, y, z):
        xn = y*z
        yn = x-y
        zn = 1-x*y
        return xn, yn, zn

class SprottLinzCAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzc_attractor_add"
    bl_label = "SprottLinzC"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    param = []
    x : get_prop("SprottLinzC", "x", 1)
    y : get_prop("SprottLinzC", "y", 0.0)
    z : get_prop("SprottLinzC", "z", 1)
    dt : get_prop("SprottLinzC", "dt", 0.02)

    def iterate(self, x, y, z):
        xn = y*z
        yn = x-y
        zn = 1-pow(x,2)
        return xn, yn, zn

class SprottLinzFAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzf_attractor_add"
    bl_label = "SprottLinzF"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzF", "a", 0.5)
    x : get_prop("SprottLinzF", "x", 0.1)
    y : get_prop("SprottLinzF", "y", 0.0)
    z : get_prop("SprottLinzF", "z", 0.0)
    dt : get_prop("SprottLinzF", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = y+z
        yn = -x+self.a*y
        zn = pow(x,2)-z
        return xn, yn, zn

class SprottLinzGAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzg_attractor_add"
    bl_label = "SprottLinzG"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzG", "a", 0.4)
    x : get_prop("SprottLinzG", "x", 1)
    y : get_prop("SprottLinzG", "y", 0.0)
    z : get_prop("SprottLinzG", "z", 0.0)
    dt : get_prop("SprottLinzG", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = self.a*x+z
        yn = x*z-y
        zn = -x+y
        return xn, yn, zn

class SprottLinzHAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzh_attractor_add"
    bl_label = "SprottLinzH"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzH", "a", 0.5)
    x : get_prop("SprottLinzH", "x", 1)
    y : get_prop("SprottLinzH", "y", 0.0)
    z : get_prop("SprottLinzH", "z", 0.0)
    dt : get_prop("SprottLinzH", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = -y+pow(z,2)
        yn = x+self.a*y
        zn = x-z
        return xn, yn, zn

class SprottLinzIAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzi_attractor_add"
    bl_label = "SprottLinzI"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzI", "a", -0.2)
    x : get_prop("SprottLinzI", "x", 0.1)
    y : get_prop("SprottLinzI", "y", 0.1)
    z : get_prop("SprottLinzI", "z", 0.1)
    dt : get_prop("SprottLinzI", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = self.a*y
        yn = x+z
        zn = x+pow(y,2)-z
        return xn, yn, zn

class SprottLinzJAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzj_attractor_add"
    bl_label = "SprottLinzJ"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzJ", "a", 2)
    x : get_prop("SprottLinzJ", "x", 0.1)
    y : get_prop("SprottLinzJ", "y", 0.1)
    z : get_prop("SprottLinzJ", "z", 0.1)
    dt : get_prop("SprottLinzJ", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = self.a*z
        yn = -self.a*y+z
        zn = -x+y+pow(y,2)
        return xn, yn, zn

class SprottLinzKAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzk_attractor_add"
    bl_label = "SprottLinzK"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzK", "a", 0.3)
    x : get_prop("SprottLinzK", "x", 0.1)
    y : get_prop("SprottLinzK", "y", 0.0)
    z : get_prop("SprottLinzK", "z", 0.0)
    dt : get_prop("SprottLinzK", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = x*y-z
        yn = x-y
        zn = x+self.a*z
        return xn, yn, zn

class SprottLinzLAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzl_attractor_add"
    bl_label = "SprottLinzL"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b"]
    a : get_prop("SprottLinzL", "a", 3.9)
    b : get_prop("SprottLinzL", "b", 0.9)
    x : get_prop("SprottLinzL", "x", 1)
    y : get_prop("SprottLinzL", "y", 1)
    z : get_prop("SprottLinzL", "z", 0.0)
    dt : get_prop("SprottLinzL", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = y+self.a*z
        yn = self.b*pow(x,2)-y
        zn = 1-x
        return xn, yn, zn

class SprottLinzMAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzm_attractor_add"
    bl_label = "SprottLinzM"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzM", "a", 1.7)
    x : get_prop("SprottLinzM", "x", 0.1)
    y : get_prop("SprottLinzM", "y", 0.1)
    z : get_prop("SprottLinzM", "z", 0.1)
    dt : get_prop("SprottLinzM", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = -z
        yn = -pow(x,2)-y
        zn = self.a+self.a*x+y
        return xn, yn, zn

class SprottLinzOAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzo_attractor_add"
    bl_label = "SprottLinzO"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzO", "a", 2.7)
    x : get_prop("SprottLinzO", "x", 0.1)
    y : get_prop("SprottLinzO", "y", 0.1)
    z : get_prop("SprottLinzO", "z", 0.1)
    dt : get_prop("SprottLinzO", "dt", 0.008)

    def iterate(self, x, y, z):
        xn = y
        yn = x-z
        zn = x+x*z+self.a*y
        return xn, yn, zn

class SprottLinzPAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzp_attractor_add"
    bl_label = "SprottLinzP"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzP", "a", 2.7)
    x : get_prop("SprottLinzP", "x", 0.1)
    y : get_prop("SprottLinzP", "y", 0.0)
    z : get_prop("SprottLinzP", "z", 0.0)
    dt : get_prop("SprottLinzP", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = self.a*y+z
        yn = -x+pow(y,2)
        zn = x+y
        return xn, yn, zn

class SprottLinzQAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzq_attractor_add"
    bl_label = "SprottLinzQ"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b"]
    a : get_prop("SprottLinzQ", "a", 3.1)
    b : get_prop("SprottLinzQ", "b", 0.5)
    x : get_prop("SprottLinzQ", "x", 0.1)
    y : get_prop("SprottLinzQ", "y", 0.1)
    z : get_prop("SprottLinzQ", "z", 0.1)
    dt : get_prop("SprottLinzQ", "dt", 0.004)

    def iterate(self, x, y, z):
        xn = -z
        yn = x-y
        zn = self.a*x+pow(y,2)+self.b*z
        return xn, yn, zn

class SprottLinzSAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottlinzs_attractor_add"
    bl_label = "SprottLinzS"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottLinzS", "a", 4)
    x : get_prop("SprottLinzS", "x", 0.1)
    y : get_prop("SprottLinzS", "y", 0.1)
    z : get_prop("SprottLinzS", "z", 0.1)
    dt : get_prop("SprottLinzS", "dt", 0.004)

    def iterate(self, x, y, z):
        xn = -x-self.a*y
        yn = x+pow(z,2)
        zn = 1+x
        return xn, yn, zn

class SprottNAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottn_attractor_add"
    bl_label = "SprottN"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a"]
    a : get_prop("SprottN", "a", 2)
    x : get_prop("SprottN", "x", 0.1)
    y : get_prop("SprottN", "y", 0.0)
    z : get_prop("SprottN", "z", 0.0)
    dt : get_prop("SprottN", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = -self.a*y
        yn = x+pow(z,2)
        zn = 1+y-self.a*z
        return xn, yn, zn

class SprottRAttractor(Object_OT_Attractor):
    bl_idname = "curve.sprottr_attractor_add"
    bl_label = "SprottR"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b"]
    a : get_prop("SprottR", "a", 0.9)
    b : get_prop("SprottR", "b", 0.4)
    x : get_prop("SprottR", "x", 0.1)
    y : get_prop("SprottR", "y", 0.0)
    z : get_prop("SprottR", "z", 0.0)
    dt : get_prop("SprottR", "dt", 0.01)

    def iterate(self, x, y, z):
        xn = self.a-y
        yn = self.b+z
        zn = x*y-z
        return xn, yn, zn

class StrizhakKawczynskiAttractor(Object_OT_Attractor):
    bl_idname = "curve.strizhakkawczynski_attractor_add"
    bl_label = "StrizhakKawczynski"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "b1", "b2", "q", "r", "ax1", "ax2", "ax3"]
    a : get_prop("StrizhakKawczynski", "a", 150)
    b : get_prop("StrizhakKawczynski", "b", 436.6)
    b1 : get_prop("StrizhakKawczynski", "b1", 3.714)
    b2 : get_prop("StrizhakKawczynski", "b2", 21.7)
    q : get_prop("StrizhakKawczynski", "q", 0.07)
    r : get_prop("StrizhakKawczynski", "r", 0.101115)
    ax1 : get_prop("StrizhakKawczynski", "ax1", 10)
    ax2 : get_prop("StrizhakKawczynski", "ax2", 11)
    ax3 : get_prop("StrizhakKawczynski", "ax3", 20)
    x : get_prop("StrizhakKawczynski", "x", 0.1)
    y : get_prop("StrizhakKawczynski", "y", 0.0)
    z : get_prop("StrizhakKawczynski", "z", 0.0)
    dt : get_prop("StrizhakKawczynski", "dt", 0.08)

    def iterate(self, x, y, z):
        xn = self.r*(y-(x-self.ax1)*(x-self.ax2)*(x-self.ax3)-self.a)
        yn = self.b-self.b1*z-self.b2*x-y
        zn = self.q*(x-z)
        return xn, yn, zn

class ThomasAttractor(Object_OT_Attractor):
    bl_idname = "curve.thomas_attractor_add"
    bl_label = "Thomas"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["b"]
    b : get_prop("Thomas", "b", 0.19)
    x : get_prop("Thomas", "x", 0.1)
    y : get_prop("Thomas", "y", 0.0)
    z : get_prop("Thomas", "z", 0.0)
    dt : get_prop("Thomas", "dt", 0.05)

    def iterate(self, x, y, z):
        xn = -self.b*x+math.sin(y)
        yn = -self.b*y+math.sin(z)
        zn = -self.b*z+math.sin(x)
        return xn, yn, zn

# TSUCS1 Attractor added by 7bitretro/googyflip 12/09/22
class TSUCS1(Object_OT_Attractor):
    bl_idname = "curve.tsucs1_attractor_add"
    bl_label = "TSUCS1"
    bl_options = {"REGISTER", "UNDO"}

    npoints : get_npoints()
    params = ["a", "b", "d", "e", "l"]
    a : get_prop("TSUCS1", "a", 40)
    b : get_prop("TSUCS1", "b", 0.833)
    d : get_prop("TSUCS1", "d", 0.5)
    e : get_prop("TSUCS1", "e", 0.65)
    l : get_prop("TSUCS1", "l", 20)
    x : get_prop("TSUCS1", "x", 0.1)
    y : get_prop("TSUCS1", "y", 0.0)
    z : get_prop("TSUCS1", "z", 0.0)
    dt : get_prop("TSUCS1", "dt", 0.001)

    def iterate(self, x, y, z):
        xn = self.a * (y - x) + self.d * x * z
        yn = self.l * y - x * z
        zn = self.b * z + x * y - self.e * x * x
        return xn, yn, zn


class MENU_MT_attractors(bpy.types.Menu):
    bl_idname= "MENU_MT_attractors"
    bl_label= "Strange Attractors"

    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_REGION_WIN"

        layout.operator("curve.act_attractor_add",text="Act")
        layout.operator("curve.aizawa_attractor_add",text="Aizawa")
        layout.operator("curve.anishenko_astakhov_attractor_add", text="Anishenko-Astakhov")
        layout.operator("curve.arneodo_attractor_add",text="Arneodo")
        layout.operator("curve.bouali_attractor_add", text="Bouali")
        layout.operator("curve.burkeshaw_attractor_add",text="Burke-Shaw")
        layout.operator("curve.chen_attractor_add",text="Chen")
        layout.operator("curve.chencelikovsky_attractor_add",text="Chen-Celikovsky")
        layout.operator("curve.chenlee_attractor_add",text="ChenLee")
        layout.operator("curve.chua_attractor_add",text="Chua")
        layout.operator("curve.chuamultii_attractor_add",text="Chua Multi I")
        layout.operator("curve.chuamultiii_attractor_add",text="Chua Multi II")
        layout.operator("curve.coullet_attractor_add", text="Coullet")
        layout.operator("curve.dadras_attractor_add", text="Dadras")
        layout.operator("curve.dequanli_attractor_add", text="Dequan-Li")
        layout.operator("curve.finance_attractor_add", text="Finance")
        layout.operator("curve.hadley_attractor_add",text="Hadley")
        layout.operator("curve.halvorsen_attractor_add",text="Halvorsen")
        layout.operator("curve.linzsprott_attractor_add",text="Sprott-Linz")
        layout.operator("curve.lorenz_attractor_add",text="Lorenz")
        layout.operator("curve.lorenzmod1_attractor_add",text="Lorenz Mod1")
        layout.operator("curve.lorenzmod2_attractor_add",text="Lorenz Mod2")
        layout.operator("curve.lotka_volterra_attractor_add",text="Lotka-Volterra")
        layout.operator("curve.moorespiegel_attractor_add",text="Moore Spiegel")
        layout.operator("curve.newjerk_attractor_add", text="New-Jerk")
        layout.operator("curve.newtonleipnik_attractor_add",text="Newton-Leipnik")
        layout.operator("curve.nosehoover_attractor_add",text="Nose-Hoover")
        layout.operator("curve.rabinovichfabrikant_attractor_add",text="Rabinovich-Fabrikant")
        layout.operator("curve.rayleighbenard_attractor_add",text="Rayleigh-Benard")
        layout.operator("curve.rikitake_attractor_add",text="Rikitake")
        layout.operator("curve.roessler_attractor_add",text="Roessler")
        layout.operator("curve.rucklidge_attractor_add",text="Rucklidge")
        layout.operator("curve.sakarya_attractor_add",text="Sakarya")
        layout.operator("curve.shimizumorioka_attractor_add",text="Shimizu-Morioka")
        layout.operator("curve.sprottd_attractor_add",text="Sprott D")
        layout.operator("curve.sprotte_attractor_add",text="Sprott E")
        layout.operator("curve.sprottlinzb_attractor_add",text="Sprott-Linz B")
        layout.operator("curve.sprottlinzc_attractor_add",text="Sprott-Linz C")
        layout.operator("curve.sprottlinzf_attractor_add",text="Sprott-Linz F")
        layout.operator("curve.sprottlinzg_attractor_add",text="Sprott-Linz G")
        layout.operator("curve.sprottlinzh_attractor_add",text="Sprott-Linz H")
        layout.operator("curve.sprottlinzi_attractor_add",text="Sprott-Linz I")
        layout.operator("curve.sprottlinzj_attractor_add",text="Sprott-Linz J")
        layout.operator("curve.sprottlinzk_attractor_add",text="Sprott-Linz K")
        layout.operator("curve.sprottlinzl_attractor_add",text="Sprott-Linz L")
        layout.operator("curve.sprottlinzm_attractor_add",text="Sprott-Linz M")
        layout.operator("curve.sprottlinzo_attractor_add",text="Sprott-Linz O")
        layout.operator("curve.sprottlinzp_attractor_add",text="Sprott-Linz P")
        layout.operator("curve.sprottlinzq_attractor_add",text="Sprott-Linz Q")
        layout.operator("curve.sprottlinzs_attractor_add",text="Sprott-Linz S")
        layout.operator("curve.sprottn_attractor_add",text="Sprott N")
        layout.operator("curve.sprottr_attractor_add",text="Sprott R")
        layout.operator("curve.strizhakkawczynski_attractor_add",text="Strizhak-Kawczynski")
        layout.operator("curve.thomas_attractor_add",text="Thomas")
        layout.operator("curve.threecellscnn_attractor_add",text="Three-Cells-CNN")
        layout.operator("curve.threelayer_attractor_add",text="3-Layer")
        layout.operator("curve.wangsun_attractor_add", text="Wang-Sun")
        layout.operator("curve.yuwang_attractor_add", text="Yu-Wang")
        layout.operator("curve.tsucs1_attractor_add", text="TSUCS1")

class INFO_MT_curve_attractor_add(bpy.types.Menu):
    bl_idname = "INFO_MT_curve_attractor_add"
    bl_label = "Attractors"

    def draw(self, context):
        layout = self.layout
        layout.operator_context = "INVOKE_REGION_WIN"

        layout.operator("curve.act_attractor_add",text="Act")
        layout.operator("curve.aizawa_attractor_add",text="Aizawa")
        layout.operator("curve.anishenko_astakhov_attractor_add", text="Anishenko-Astakhov")
        layout.operator("curve.arneodo_attractor_add",text="Arneodo")
        layout.operator("curve.bouali_attractor_add", text="Bouali")
        layout.operator("curve.burkeshaw_attractor_add",text="Burke-Shaw")
        layout.operator("curve.chen_attractor_add",text="Chen")
        layout.operator("curve.chencelikovsky_attractor_add",text="Chen-Celikovsky")
        layout.operator("curve.chenlee_attractor_add",text="ChenLee")
        layout.operator("curve.chua_attractor_add",text="Chua")
        layout.operator("curve.chuamultii_attractor_add",text="Chua Multi I")
        layout.operator("curve.chuamultiii_attractor_add",text="Chua Multi II")
        layout.operator("curve.coullet_attractor_add", text="Coullet")
        layout.operator("curve.dadras_attractor_add", text="Dadras")
        layout.operator("curve.dequanli_attractor_add", text="Dequan-Li")
        layout.operator("curve.finance_attractor_add", text="Finance")
        layout.operator("curve.hadley_attractor_add",text="Hadley")
        layout.operator("curve.halvorsen_attractor_add",text="Halvorsen")
        layout.operator("curve.linzsprott_attractor_add",text="Sprott-Linz")
        layout.operator("curve.lorenz_attractor_add",text="Lorenz")
        layout.operator("curve.lorenzmod1_attractor_add",text="Lorenz Mod1")
        layout.operator("curve.lorenzmod2_attractor_add",text="Lorenz Mod2")
        layout.operator("curve.lotka_volterra_attractor_add",text="Lotka-Volterra")
        layout.operator("curve.moorespiegel_attractor_add",text="Moore Spiegel")
        layout.operator("curve.newjerk_attractor_add", text="New-Jerk")
        layout.operator("curve.newtonleipnik_attractor_add",text="Newton-Leipnik")
        layout.operator("curve.nosehoover_attractor_add",text="Nose-Hoover")
        layout.operator("curve.rabinovichfabrikant_attractor_add",text="Rabinovich-Fabrikant")
        layout.operator("curve.rayleighbenard_attractor_add",text="Rayleigh-Benard")
        layout.operator("curve.rikitake_attractor_add",text="Rikitake")
        layout.operator("curve.roessler_attractor_add",text="Roessler")
        layout.operator("curve.rucklidge_attractor_add",text="Rucklidge")
        layout.operator("curve.sakarya_attractor_add",text="Sakarya")
        layout.operator("curve.shimizumorioka_attractor_add",text="Shimizu-Morioka")
        layout.operator("curve.sprottd_attractor_add",text="Sprott D")
        layout.operator("curve.sprotte_attractor_add",text="Sprott E")
        layout.operator("curve.sprottlinzb_attractor_add",text="Sprott-Linz B")
        layout.operator("curve.sprottlinzc_attractor_add",text="Sprott-Linz C")
        layout.operator("curve.sprottlinzf_attractor_add",text="Sprott-Linz F")
        layout.operator("curve.sprottlinzg_attractor_add",text="Sprott-Linz G")
        layout.operator("curve.sprottlinzh_attractor_add",text="Sprott-Linz H")
        layout.operator("curve.sprottlinzi_attractor_add",text="Sprott-Linz I")
        layout.operator("curve.sprottlinzj_attractor_add",text="Sprott-Linz J")
        layout.operator("curve.sprottlinzk_attractor_add",text="Sprott-Linz K")
        layout.operator("curve.sprottlinzl_attractor_add",text="Sprott-Linz L")
        layout.operator("curve.sprottlinzm_attractor_add",text="Sprott-Linz M")
        layout.operator("curve.sprottlinzo_attractor_add",text="Sprott-Linz O")
        layout.operator("curve.sprottlinzp_attractor_add",text="Sprott-Linz P")
        layout.operator("curve.sprottlinzq_attractor_add",text="Sprott-Linz Q")
        layout.operator("curve.sprottlinzs_attractor_add",text="Sprott-Linz S")
        layout.operator("curve.sprottn_attractor_add",text="Sprott N")
        layout.operator("curve.sprottr_attractor_add",text="Sprott R")
        layout.operator("curve.strizhakkawczynski_attractor_add",text="Strizhak-Kawczynski")
        layout.operator("curve.thomas_attractor_add",text="Thomas")
        layout.operator("curve.threecellscnn_attractor_add",text="Three-Cells-CNN")
        layout.operator("curve.threelayer_attractor_add",text="3-Layer")
        layout.operator("curve.wangsun_attractor_add", text="Wang-Sun")
        layout.operator("curve.yuwang_attractor_add", text="Yu-Wang")
        layout.operator("curve.tsucs1_attractor_add", text="TSUCS1")

def menu_func(self, context):
    self.layout.menu("INFO_MT_curve_attractor_add", icon="PLUGIN")

def menu_func(self, context):
    self.layout.menu("MENU_MT_attractors", icon="PLUGIN")

#needed as workaround
bpy.types.VIEW3D_MT_curve_add.append(menu_func)

#All classes used
classes = (
            Object_OT_Attractor,
            MENU_MT_attractors,
#            LorenzAttractor,
#            CoulletAttractor,
#            ActAttractor
CoulletAttractor,
LorenzAttractor,
RoesslerAttractor,
AizawaAttractor,
ActAttractor,
ThreeCellsCNNAttractor,
ArneodoAttractor,
AnishenkoAstakhovAttractor,
BoualiAttractor,
BurkeShawAttractor,
ChenAttractor,
LotkaVolterraAttractor,
MooreSpiegelAttractor,
RikitakeAttractor,
RabinovichFabrikantAttractor,
ThreeLayerAttractor,
ChuaAttractor,
ChuaMultiIAttractor,
ChuaMultiIIAttractor,
ChenLeeAttractor,
FinanceAttractor,
ChenCelikovskyAttractor,
DadrasAttractor,
DequanLiAttractor,
YuWangAttractor,
HadleyAttractor,
HalvorsenAttractor,
LinzSprottAttractor,
LorenzMod1Attractor,
LorenzMod2Attractor,
NewtonLeipnikAttractor,
NoseHooverAttractor,
RayleighBenardAttractor,
RucklidgeAttractor,
WangSunAttractor,
NewJerkAttractor,
SakaryaAttractor,
Duffing,
ShimizuMoriokaAttractor,
SprottDAttractor,
SprottEAttractor,
SprottLinzBAttractor,
SprottLinzCAttractor,
SprottLinzFAttractor,
SprottLinzGAttractor,
SprottLinzHAttractor,
SprottLinzIAttractor,
SprottLinzJAttractor,
SprottLinzKAttractor,
SprottLinzLAttractor,
SprottLinzMAttractor,
SprottLinzOAttractor,
SprottLinzPAttractor,
SprottLinzQAttractor,
SprottLinzSAttractor,
SprottNAttractor,
SprottRAttractor,
StrizhakKawczynskiAttractor,
ThomasAttractor,
TSUCS1
)

#Register and unregister all classes
register, unregister = bpy.utils.register_classes_factory(classes)

# This allows you to run the script directly from Blender's Text editor
# to test the add-on without having to install it.
if __name__ == "__main__":
  register()
