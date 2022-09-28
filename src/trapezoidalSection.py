from argparse import ArgumentError
from src.perf import Performance
from src.utils import G
from src.irregularSection import IrregularSection
from scipy.optimize import newton

import numpy as np

class TrapezoidalSection(IrregularSection):

    def __init__(self, x, z, b0, b0_min, s_right=None, s_left=None, f_right=None, f_left=None, z_min=None, y_max=None , up_section=None, down_section=None, granulometry=None, manning=None, K_over_tauc=None, tauc_over_rho=None):
        self.__b0 = b0
        if b0_min > b0:
            raise(ArgumentError("b0_min cannot be greater than b0 for a trapezoidal section"))
        self.__b0_min = b0_min
        if s_right != None:
            self.__f_right = 1/s_right
        elif f_right != None:
            self.__f_right = f_right
        else:
            raise(ArgumentError("you must give slope or 'fruit' (reciprocal of slope) when you define a trapezoidal section."))
        if s_left != None:
            self.__f_left = 1/s_left
        elif f_left != None:
            self.__f_left = f_left
        else:
            raise(ArgumentError("you must give slope or 'fruit' (reciprocal of slope) when you define a trapezoidal section."))
        
        self.__h = 0
        self.__y_max = 1000 if y_max==None or y_max <= 0 else y_max # MAX_INT
        # points = [(0, self.__y_max), (self.__f_right*self.__y_max,0), (self.__f_right*self.__y_max+b0, 0), (self.__f_right*self.__y_max+self.__f_left*self.__y_max+b0, self.__y_max)]
        points = self.get_points(self.__b0, self.__h)
        super().__init__(points, x, z, z_min=z_min, up_section=up_section, down_section=down_section, granulometry=granulometry, manning=manning, K_over_tauc=K_over_tauc, tauc_over_rho=tauc_over_rho)
        
    def interp_as_up_section(self, other_section, x=None):
        interpolated_section = super().interp(self, other_section, x=x)
        if x==None:
            x = 0.5*(self.get_x()+other_section.get_x())
        interpolated_section.__b0 = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_b0(), other_section.get_b0()])
        interpolated_section.__b0_min = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_b0_min(), other_section.get_b0_min()])
        interpolated_section.__f_right = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_f_right(), other_section.get_f_right()])
        interpolated_section.__f_left = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_f_left(), other_section.get_f_left()])
        return interpolated_section

    def interp_as_down_section(self, other_section, x=None):
        interpolated_section = super().interp(other_section, self, x=x)
        if x==None:
            x = 0.5*(self.get_x()+other_section.get_x())
        interpolated_section.__b0 = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_b0(), other_section.get_b0()])
        interpolated_section.__b0_min = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_b0_min(), other_section.get_b0_min()])
        interpolated_section.__f_right = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_f_right(), other_section.get_f_right()])
        interpolated_section.__f_left = np.interp(x, [self.get_x(), other_section.get_x()], [self.get_f_left(), other_section.get_f_left()])
        return interpolated_section

    def copy(self):
        """return a safe copy of this section"""
        return TrapezoidalSection(x=self.get_x(), z=self.get_z(), b0=self.get_b0(), b0_min=self.get_b0_min(), f_right=self.get_f_right(), f_left=self.get_f_left(), z_min=self.get_z_min(), y_max=self.get_y_max() , up_section=self.get_up_section(), down_section=self.get_down_section(), granulometry=self.get_granulometry(), manning=self.get_manning(), K_over_tauc=None, tauc_over_rho=None)

    def update_bottom(self, Q, y_up, y, y_down, Qs_in, dt, law):
        z_min = self.get_z_min()
        z = self.get_z()
        inter_section_up = self.interp_as_down_section(self.get_up_section())
        inter_section_down = self.interp_as_up_section(self.get_down_section())
        b0_up, h_up = inter_section_up.get_b0(), inter_section_up.get_h()
        b0_down, h_down = inter_section_down.get_b0(), inter_section_down.get_h()
        b0, h = self.get_b0(), self.get_h()
        dx_up = abs(inter_section_up.get_x()-self.get_x())
        dx_down = abs(inter_section_down.get_x()-self.get_x())
        # Numerical convention for upstream and downstream section, not very important.
        if dx_up==0:
            dx_up = dx_down
        if dx_down==0:
            dx_down = dx_up
        ######
        def func_to_zero(z_new, vol):
            delta_z =  z_new - z
            b0_new_up, h_new_up = inter_section_up.get_new_b0_and_h(delta_z)
            b0_new_down, h_new_down = inter_section_down.get_new_b0_and_h(delta_z)
            b0_new, h_new = self.get_new_b0_and_h(delta_z)

            delta_h_up = h_new_up - h_up
            delta_h_down = h_new_down - h_down
            delta_h = h_new - h
            V_up_1 = 0.5*dx_up*(-delta_h*b0 - delta_h_up*b0_up)
            V_down_1 = 0.5*dx_down*(-delta_h*b0 - delta_h_down*b0_down)
            V_up_2 = 0.25*dx_up*((delta_z+delta_h)*(b0 + b0_new)+(delta_z+delta_h_up)*(b0_up+b0_new_up))
            V_down_2 = 0.25*dx_down*((delta_z+delta_h)*(b0 + b0_new)+(delta_z+delta_h_down)*(b0_down+b0_new_down))
            return V_up_1+V_up_2+V_down_1+V_down_2 - vol

        Qs_out = law.compute_Qs(self, Q, y, y_down)
        Qs_out_max = -func_to_zero(z_min, 0)/dt
        if self.is_downstream():
            Qs_out = Qs_in
        elif Qs_out > Qs_out_max:
            Qs_out = Qs_out_max
        
        vol = (Qs_in-Qs_out)*dt
        z_new = newton(func_to_zero, x0=z, args=[vol])
        
        b0_new, h_new = self.get_new_b0_and_h(z_new-z)
        self.set_z(z_new)
        self.set_b0(b0_new)
        self.set_h(h_new)
        self.__points = self.get_points(b0_new, h_new)

        return Qs_out

    def get_stored_volume(self):
        z_min = self.get_z_min()
        z = self.get_z()
        inter_section_up = self.interp_as_down_section(self.get_up_section())
        inter_section_down = self.interp_as_up_section(self.get_down_section())
        b0_up, h_up = inter_section_up.get_b0(), inter_section_up.get_h()
        b0_down, h_down = inter_section_down.get_b0(), inter_section_down.get_h()
        b0, h = self.get_b0(), self.get_h()
        dx_up = abs(inter_section_up.get_x()-self.get_x())
        dx_down = abs(inter_section_down.get_x()-self.get_x())
        # Numerical convention for upstream and downstream section, not very important.
        if dx_up==0:
            dx_up = dx_down
        if dx_down==0:
            dx_down = dx_up
        ######
        delta_z = z-z_min
        b0_new_up, h_new_up = inter_section_up.get_new_b0_and_h(delta_z)
        b0_new_down, h_new_down = inter_section_down.get_new_b0_and_h(delta_z)
        b0_new, h_new = self.get_new_b0_and_h(delta_z)

        delta_h_up = h_new_up - h_up
        delta_h_down = h_new_down - h_down
        delta_h = h_new - h
        V_up_1 = 0.5*dx_up*(delta_h*b0 + delta_h_up*b0_up)
        V_down_1 = 0.5*dx_down*(delta_h*b0 + delta_h_down*b0_down)
        V_up_2 = 0.25*dx_up*((delta_z-delta_h)*(b0 + b0_new)+(delta_z-delta_h_up)*(b0_up+b0_new_up))
        V_down_2 = 0.25*dx_down*((delta_z-delta_h)*(b0 + b0_new)+(delta_z-delta_h_down)*(b0_down+b0_new_down))
        return V_up_1+V_up_2+V_down_1+V_down_2

    # getters and setters

    def get_b0(self):
        return self.__b0

    def set_b0(self, b0):
        if b0 < self.get_b0_min():
            raise ValueError("in a trapezoidal section b0 must be >= b0_min")
        self.__b0 == b0

    def get_b0_min(self):
        return self.__b0_min

    def get_f_right(self):
        return self.__f_right

    def get_f_left(self):
        return self.__f_left

    def get_h(self):
        """
        \     /
         \   /
          |_| = h
        """
        return self.__h

    def set_h(self, h):
        if h < 0:
            raise ValueError("in a trapezoidal section h muste be >= 0")
        self.__h = h 

    def get_wet_section(self, y):
        points = self.get_points()
        h = self.get_h()    
        if y >= self.get_y_max():
            print("WARNING : Water depth has gone upper than the maximum one")
            return points
        elif h == 0:
            return [(self.get_f_right()*(self.__y_max - y), y)] + points[1:-1] + [(self.get_f_right()*self.__y_max+self.get_b0()+y*self.get_f_left(), y)]
        else:
            if y < h:
                return [(points[2][0], y)] + points[2:-2] + [(points[-2][0], y)]
            else:
                return [(self.get_f_right()*(self.__y_max - (y-h)), y-h)] + points[1:-1] + [(self.get_f_right()*self.__y_max+self.get_b0()+(y-h)*self.get_f_left(), y-h)]
                

    def get_b(self, y=0, wet_points=None):
        return self.get_b0() + y*self.get_f_right() + y*self.get_f_left()

    def set_b0(self, b0):
        if b0 <= 0:
            raise ValueError("width b0 can not be lower than 0")
        self.__b0 = b0

    def get_S(self, y, wet_points=None):
        return (self.get_b0()+self.get_b(y, wet_points))*0.5*y
    
    def get_P(self, y, wet_points=None):
        return self.get_b0() + y*np.sqrt(1+self.get_f_right()**2) + y*np.sqrt(1+self.get_f_left()**2)

    def get_R(self, y, wet_points=None):
        return self.get_S(y) / self.get_P(y)

    def get_V(self, Q, y, wet_points=None):
        return Q / self.get_S(y)

    def get_H(self, Q, y, wet_points=None):
        return self.get_z() + y + (self.get_V(Q, y)**2)/(2*G)

    def get_Hs(self, Q, y, wet_points=None):
        return y + (self.get_V(Q, y)**2)/(2*G)

    def get_Fr(self, Q, y, wet_points=None):
        return self.get_V(Q, y)/((G * y)**0.5)

    def __get_A_for_coussot(self, y):
        return 1.93 - 0.6 * np.arctan((0.4*y/self.get_b0())**20)

    @Performance.measure_perf
    def get_yc(self, Q):
        coefs = np.array([(self.get_f_right()+self.get_f_left())**2, 4*self.get_b0()*(self.get_f_right()+self.get_f_left()), 4*self.get_b0()**2, 0, 0, -4*Q**2/G])
        r = np.roots(coefs)
        real = r.real[abs(r.imag)<1e-5]
        real_positive = real[real>0]
        if real_positive.size == 0:
            super().get_yc(Q)
        else:
            return min(real_positive)


    # Operators overloading

    def __str__(self):
        return f'Trapezoidal section : x={self.get_x()}, z={self.get_z()}, b={self.__b}, f_left={self.get_f_left()}, f_right={self.get_f_right()}'

    # EXPERIENCES

    def get_new_b0_and_h(self, delta_z):
        if self.get_z() + delta_z < self.get_z_min():
            delta_z = self.get_z_min() - self.get_z()
        h = self.get_h()
        b0 = self.get_b0()
        b0_min = self.get_b0_min()
        if h > 0:
            if delta_z < h:
                return b0, h-delta_z
            else:
                delta_z_1 = delta_z-h
                new_b0 = b0 + delta_z_1*(self.get_f_right() + self.get_f_left())
                return new_b0, 0
        else:
            new_b0 = b0 + delta_z*(self.get_f_right() + self.get_f_left())
            if new_b0 >= b0_min:
                return new_b0, 0
            else:
                delta_z_1 = (b0_min-b0)/(self.get_f_right() + self.get_f_left()) if (self.get_f_right() + self.get_f_left() != 0) else 0
                return b0_min, abs(delta_z-delta_z_1)

    def get_points(self, new_b0=None, new_h=None):
        new_b0 = new_b0 if new_b0 != None else self.get_b0()
        new_h = new_h if new_h != None else self.get_h()
        if new_h==0:
            return [(0, self.__y_max), (self.__f_right*self.__y_max,0), (self.__f_right*self.__y_max+new_b0, 0), (self.__f_right*self.__y_max+self.__f_left*self.__y_max+new_b0, self.__y_max)]
        else:
            return [(0, self.__y_max), (self.__f_right*self.__y_max, new_h), (self.__f_right*self.__y_max,0), (self.__f_right*self.__y_max+new_b0, 0), (self.__f_right*self.__y_max+new_b0, new_h), (self.__f_right*self.__y_max+self.__f_left*self.__y_max+new_b0, self.__y_max)]

