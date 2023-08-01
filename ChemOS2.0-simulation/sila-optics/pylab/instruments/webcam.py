import pyvisa.constants as pv_const
import time


class RTSP(object):
    def __init__(self, visa, location, *args, **kwargs):
        self.cam = cv2.VideoCapture(visa)
        while not self.cam.isOpened():
            time.sleep(1)

        self.location = location
        self.crop     = kwargs.get('crop', None)
        self.quality  = kwargs.get('quality', 100)
        self.shrink   = kwargs.get('shrink', 1.0)

    def __del__(self):
        self.manager.uv_led.off()
        del self.manager

    # ---- override methods ----
    def write(self, command):
        return

    def ask(self, command):
        return self.read()

    def read(self):
        ret, frame = self.cam.read()
        return frame

    
    # ---- heating ----
    def save_frame(self, fname, process = True):
        if process:
            frame = self.get_frame()
        else:
            frame = self.read()
        full_path = os.path.join(self.location, fname)
        cv2.imwrite(full_path, frame, [int(cv2.IMWRITE_JPEG_QUALITY), self.quality])

    def get_frame(self):
        frame = self.read()

        if self.crop:
            y1, y2, x1, x2 = self.crop
            frame = frame[y1:y2, x1:x2]
        frame = cv2.resize(frame, None, fx=shrink, fy=shrink)




