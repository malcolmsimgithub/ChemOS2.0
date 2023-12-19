import numpy
import cv2
from tkinter import filedialog
import math
import os

class ImageHandler():
    
    def __init__(self):
        self.filetypes =   (('PNG files', '*.png'),
                            ('JPEG files', '*.jpeg, *.jpg, *.jpe'),
                            ('TIFF files', '*.tiff, *.tif'),
                            ('bitmaps files', '*.bmp, *.dib'),
                            ('JPEG 2000 files', '*.jp2'),
                            ('WebP', '*.webp'),
                            ('Portable image format', '*.pbm, *.pgm, *.ppm *.pxm, *.pnm'),
                            ('PFM files', '*.pfm'),
                            ('Sun rasters', '*.sr, *.ras'),
                            ('OpenEXR Image files', '*.exr'),
                            ('Radiance HDR', '*.hdr, *.pic'))


    def open_images(self, filelist, grayscale = False):
        im_list = []
        for filename in filelist:
            im = cv2.imread(filename)
            if grayscale:
                im= cv2.cvtColor(im, cv2.COLOR_BGR2GRAY)
            im_list.append(im)
        return im_list

    def open_images_dialog(self, grayscale = False):

        filelist = filedialog.askopenfilenames(filetypes = self.filetypes)
        print(filelist)
        im_list = self.open_images(filelist, grayscale=grayscale)
        return im_list, filelist


    def crop_image(self, img, mergin):

        if mergin[1] == 0 and mergin[3] == 0:
            return img[mergin[0]:, mergin[2]:, : ]
        elif mergin[1] == 0:
            return img[mergin[0]:, mergin[2]:-mergin[3], : ]
        elif mergin[3] == 0:
            return img[mergin[0]:-mergin[1], mergin[2]:, : ] 
        else:
            return img[mergin[0]:-mergin[1], mergin[2]:-mergin[3], : ]


    def tile_img(self, h_num, save_filename, mergin = None, grayscale = False, del_files = False, filelist = None):
        
        if filelist == None:
            im_list, filelist = self.open_images_dialog(grayscale=False)
        else:
            im_list = self.open_images(filelist, grayscale=False)

        #defines the tile
        im_num = len(im_list)
        v_num = math.ceil(im_num / h_num)

        #cut mergin from images
        if mergin:       
            im_list = [self.crop_image(img, mergin) for img in im_list] 

        #resize image to the minimum size in the list
        h_min = min(im.shape[0] for im in im_list)
        v_min = min(im.shape[1] for im in im_list)
        im_list = [cv2.resize(im, (v_min, h_min), interpolation=cv2.INTER_CUBIC) for im in im_list]

        #add white image to the image gap
        for i in range(v_num * h_num - im_num):
            im_list.append(im_list[0]*0 + 255)
        
        #tile image
        im = cv2.vconcat([cv2.hconcat([im_list[i+h_num*j] for i in range(h_num)]) for j in range(v_num)]) 

        #savefig
        cv2.imwrite(save_filename, im)

        #delete orignal images
        if del_files:
            for f in filelist:
                if f != save_filename:
                    try:
                        os.remove(f)
                    except:
                        pass    
        return im

if __name__ == '__main__':

    IM = ImageHandler()

    IM.tile_img(h_num = 4, save_filename = 'D:/temp_for_covid19/HPLCMS/20200707_splitter_test/without_splitter/MS_peaks.png', mergin = [0,0,0,0])


# im_v = cv2.hconcat(im)
# cv2.imwrite('D:/temp_for_covid19/MS_data/test.png', im_v)