import openslide
import torch
import torchvision.transforms as transforms
import numpy as np
import matplotlib.pyplot as plt
from xml.dom import minidom
from options import options
import options
import math

options = options.get_options()

#wsi_path = "WSI/TCGA-EB-A44Q-01Z-00-DX1.E7EA9878-E2B1-4768-8157-4459BDE753F0.svs"
#name = "TCGA-GN-A8LN"

def read_xml(xml_path):
    annotated_file = minidom.parse(xml_path)
    coordinates = annotated_file.getElementsByTagName("Coordinate")
    pairs = []
    for i in range(len(coordinates)):
        if i%2 == 0:
            x1 = int(float(coordinates[i].attributes["X"].value))
            y1 = int(float(coordinates[i].attributes["Y"].value))
            x2 = int(float(coordinates[i+1].attributes["X"].value))
            y2 = int(float(coordinates[i+1].attributes["Y"].value))
            pair = [(x1,y1), (x2,y2)]
            pairs.append(pair)
    return pairs

def rotation_remap(x, y, theta, img_size):
    theta = math.radians(-theta)
    x_center, y_center = img_size[0]/2, img_size[1]/2
    x, y = x-x_center, y-y_center
    x_remap = x*math.cos(theta) - y*math.sin(theta)
    y_remap = x*math.sin(theta) + y*math.cos(theta)
    x_remap, y_remap = x_remap+x_center, y_remap+y_center
    x_remap = max(x_remap, 0)
    x_remap = min(x_remap, img_size[0])
    y_remap = max(y_remap, 0)
    y_remap = min(y_remap, img_size[1])
    return x_remap, y_remap


def get_biopsy(start_location, direct_location, wsi, pixel_length, pixel_radius):
    mpp = (float(wsi.properties['openslide.mpp-x']) + float(wsi.properties['openslide.mpp-x']))/2
    x1, y1 = start_location
    x2, y2 = direct_location

    x_direction = ((x2-x1)/abs(x2-x1) if x2-x1 != 0 else 0)
    y_direction = ((y2-y1)/abs(y2-y1) if y2-y1 != 0 else 0)


    x3 = x1 + pixel_length * x_direction
    y3 = y1 + pixel_length * y_direction


    if x_direction != 0 and y_direction != 0:

        tan = (x2-x1)/(y2-y1)
        angle = (math.atan(tan) if x_direction*y_direction > 0 else -math.atan(tan))
        angle_deg = math.degrees( angle )

        x3 = x1 + pixel_length * math.sin(angle) * x_direction
        y3 = y1 + pixel_length * math.cos(angle) * y_direction

        diagonal_length = (pixel_radius**2 + pixel_length**2)**(1/2)
        center_x = x1 + 0.5 * pixel_length * math.sin(angle) * x_direction
        center_y = y1 + 0.5 * pixel_length * math.cos(angle) * y_direction
        crop_x = center_x - 0.5 * pixel_length
        crop_y = center_y - 0.5 * pixel_length
        crop_location = ((int(crop_x), int(crop_y)))
        box_size = (int(diagonal_length)+2, int(diagonal_length)+2)
        img = wsi.read_region(location=crop_location, size=box_size, level=0)
        x1_before_rotate = x1 - crop_location[0]
        y1_before_rotate = y1 - crop_location[1]
        x3_before_rotate = x3 - crop_location[0]
        y3_before_rotate = y3 - crop_location[1]
        if x_direction > 0 and y_direction > 0:
            rotation_angle = 90-angle_deg
        elif x_direction > 0 and y_direction < 0:
            rotation_angle = angle_deg-90
        elif x_direction < 0 and y_direction > 0:
            rotation_angle = 90+angle_deg
        elif x_direction < 0 and y_direction < 0:
            rotation_angle =-90-angle_deg


        x1_after_rotate, y1_after_rotate = rotation_remap(x=x1_before_rotate,
                                                         y=y1_before_rotate,
                                                         theta=rotation_angle,
                                                         img_size=box_size)
        x3_after_rotate, y3_after_rotate = rotation_remap(x=x3_before_rotate,
                                                         y=y3_before_rotate,
                                                         theta=rotation_angle,
                                                         img_size=box_size)


        rotated_img = img.rotate(rotation_angle)
        top = max(min(y1_after_rotate - 0.5*pixel_radius, y1_after_rotate + 0.5*pixel_radius, 
                   y3_after_rotate - 0.5*pixel_radius, y3_after_rotate + 0.5*pixel_radius),0)
        bottom = min(max(y1_after_rotate - 0.5*pixel_radius, y1_after_rotate + 0.5*pixel_radius, 
                   y3_after_rotate - 0.5*pixel_radius, y3_after_rotate + 0.5*pixel_radius), box_size[1])
        left = max(min(x1_after_rotate, x3_after_rotate), 0)
        right = min(max(x1_after_rotate, x3_after_rotate), box_size[0])
        sample = rotated_img.crop((left, top, right, bottom))
    else:
        if x_direction == 0 and y_direction > 0:
            sample = wsi.read_region(location=(int(x1-pixel_radius*0.5),int(y1)), size=(int(pixel_length), int(pixel_length)), level=0)
            sample = sample.rotate(90)
            sample = sample.crop((0, (pixel_length-pixel_radius)/2, pixel_length, (pixel_length+pixel_radius)/2))
        elif x_direction == 0 and y_direction < 0:
            sample = wsi.read_region(location=(int(x3-pixel_radius*0.5),int(y3)), size=(int(pixel_length), int(pixel_length)), level=0)
            sample = sample.crop((0, (pixel_length-pixel_radius)/2, pixel_length, (pixel_length+pixel_radius)/2))
        elif x_direction > 0 and y_direction == 0:
            sample = wsi.read_region(location=(int(x1),int(y1-pixel_radius*0.5)), size=(int(pixel_length), int(pixel_length)), level=0)
            sample = sample.crop((0, (pixel_length-pixel_radius)/2, pixel_length, (pixel_length+pixel_radius)/2))
        elif x_direction < 0 and y_direction == 0:
            sample = wsi.read_region(location=(int(x3),int(y3-pixel_radius*0.5)), size=(int(pixel_length), int(pixel_length)), level=0)
            sample = sample.rotate(180)
            sample = sample.crop((0, (pixel_length-pixel_radius)/2, pixel_length, (pixel_length+pixel_radius)/2))
  
            
    return sample

def make_biopsy(wsi_path, start_location, direct_location):
    wsi = openslide.OpenSlide(wsi_path)
    mpp = (float(wsi.properties['openslide.mpp-x']) + float(wsi.properties['openslide.mpp-x']))/2
    biopsy_length = np.random.normal(options['length'], options['length_sd']) * 1000 # convert to micron
    biopsy_radius = np.random.normal(options['radius'], options['radius_sd']) * 1000
    pixel_length = int(biopsy_length/mpp)
    pixel_radius = int(biopsy_radius/mpp)
    sample = get_biopsy(start_location=start_location, 
                 direct_location=direct_location, 
                 wsi=wsi, 
                 pixel_length=pixel_length, 
                 pixel_radius=pixel_radius)
    return sample

def make_slide(samples):
    n_core = len(samples)
    max_length = 0
    total_radius = 0
    for sample in samples:
        total_radius += sample.size[1]
        if sample.size[0] > max_length:
            max_length = sample.size[0]
    slide_size = (3, int(total_radius*2.0), int(max_length*1.5))
    img = torch.normal(0.9,0.01, size=slide_size) # calling background
    space = int(total_radius/n_core)
    location_y = int((slide_size[2]-max_length)/2)
    location_x = int(space/2)
    to_tensor = transforms.ToTensor()
    for sample in samples:
        img[:,location_x:(location_x+sample.size[1]), location_y:(location_y+sample.size[0])] = to_tensor(sample)[:3]
        location_x = int(location_x + sample.size[1] + space)
    to_img = transforms.ToPILImage()
    img = to_img(img)
    return img

def get_biopsy_from_wsi(wsi_path, xml_path):
    pairs = read_xml(xml_path)
    samples = []
    for pair in pairs:
        start_location, direct_location = pair
        sample = make_biopsy(wsi_path, start_location=start_location, direct_location=direct_location)
        samples.append(sample)
    slide = make_slide(samples)
    return slide
