
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

path = os.path.dirname(os.path.abspath(__file__))

#color matching function from http://cvrl.ioo.ucl.ac.uk/cmfs.htm (CIE 2006 2 Degrees)
cmf_filename = '%s/lin2012xyz2e_1_7sf.csv' %path 



class Spec_to_Color():
    def __init__(self):
        self.cmf =  np.transpose(np.loadtxt(cmf_filename, delimiter=','))

        #matrix used to convert RGB to XYZ (https://colorcalculations.wordpress.com/rgb-color-spaces/#RGBspaces)
        self.M_CIE = np.asarray([[0.488718,0.310680,0.200602],[0.176204,0.812985,0.010811],[0.000000,0.010205,0.989795]])
        self.MI_CIE =  np.linalg.inv(self.M_CIE)
        self.M_NTSC = np.asarray([[0.601404,0.174950,0.196503], [0.296214,0.591499,0.112287], [0.000000,0.066648,1.094800]])
        self.M_sRGB =  np.asarray([[0.412529,0.358164,0.177404],[0.212710,0.716328,0.070961],[0.019337,0.119388,0.934326]])
        self.M_Adobe =  np.asarray([[0.575930,0.188922,0.183244],[0.287965,0.638737,0.073298],[0.035996,0.071970,0.965085]])
        self.M_Apple =  np.asarray([[0.449941,0.316877,0.181278],[0.244768,0.673364,0.081868],[0.025197,0.141463,0.906392]])

    def interpolate(self, spec, kind = 'cubic'):
        f = interpolate.interp1d(spec[0], spec[1], kind=kind)
        
        return [self.cmf[0], f(self.cmf[0])]

    def spec_to_xyz(self, spec, interpolate='cubic'):

        i_spec = self.interpolate(spec, kind = 'cubic')

        XYZ = [np.sum(self.cmf[1] * i_spec[1]), np.sum(self.cmf[2] * i_spec[1]), np.sum(self.cmf[3] * i_spec[1])]
        den = sum(XYZ)
        xyz = XYZ/den

        return xyz

    def spec_to_rgb(self, spec, interpolate='cubic', color_space = 'sRGB'):

        xyz = self.spec_to_xyz(spec, interpolate=interpolate)
        rgb = self.xyz_to_rgb(xyz, color_space = color_space)

        return rgb

    def xyz_to_rgb(self, xyz, color_space = 'sRGB'):

        if color_space == 'CIE':
            MI = np.linalg.inv(self.M_CIE)
        elif color_space == 'NTSC':
            MI = np.linalg.inv(self.M_NTSC)
        elif color_space == 'sRGB':
            MI = np.linalg.inv(self.M_sRGB)
        elif color_space == 'Adobe':
            MI = np.linalg.inv(self.M_Adobe)
        elif color_space == 'Apple':
            MI = np.linalg.inv(self.M_Apple)

        rgb = MI.dot(xyz)

        # We're not in the RGB gamut: approximate by desaturating
        if np.any(rgb < 0):
            w = - np.min(rgb)
            rgb += w

        # Normalize the rgb vector (since we don't care about the power)
        if not np.all(rgb==0):
            rgb /= np.max(rgb)
        
    
        return self.rgb_to_hex(rgb)

    def rgb_to_hex(self, rgb):
        return (rgb*255).astype(int)

    def hex_to_str(self, rgb):
        text = '#'
        for v in rgb:
            text += format(v, '02x')
        return text

if __name__ == '__main__':

    SC = Spec_to_Color()
    spec = np.transpose(np.loadtxt('%s/201030_12-14_PL.txt' %path, usecols=[0,1,2,3]))

    spec = [spec[0], spec[3]]

    xyz = SC.spec_to_xyz(spec)

    print(xyz)

    rgb = SC.xyz_to_rgb(xyz)
    print(rgb, SC.hex_to_str(rgb))

    # colar_spaces = ['CIE', 'NTSC', 'sRGB', 'Adobe', 'Apple']

    # for space in colar_spaces:

    #     rgb = SC.xyz_to_rgb(xyz, color_space=space)

    #     print(rgb)
    #     print(SC.hex_to_str(rgb))

