import openslide
import torch
import torchvision.transforms as transforms
import numpy as np
import matplotlib.pyplot as plt
from xml.dom import minidom
from options import options
import options

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

def get_biopsy(start_location, direct_location, wsi, pixel_length, pixel_radius):
    mpp = float(wsi.properties['openslide.mpp-x'])
    points = [start_location, direct_location]

    x1, y1 = start_location
    x2, y2 = direct_location
    distance_x = (x2-x1)
    distance_y = (y2-y1)
    
    if distance_x > 0:
        ratio = distance_y/distance_x
        current_location = start_location
        # create biopy
        diagonal = 0
        current_img = (torch.zeros(3,pixel_radius,1) if abs(ratio)<=1 else torch.zeros(3,1,pixel_radius))
        size = ((1, pixel_radius) if abs(ratio)<=1 else (pixel_radius,1))
        dim = (2 if abs(ratio)<=1 else 1)
        move = (1 if y1 < y2 else -1)
        to_tensor = transforms.ToTensor()
        to_img = transforms.ToPILImage()
        while diagonal <= pixel_length:
            current_location = ((current_location[0] + 1, current_location[1] + ratio) if abs(ratio)<=1 else (current_location[0] + move/ratio, current_location[1] + move))
            each_img = wsi.read_region(location=(int(current_location[0]), int(current_location[1])),size=size, level=0)
            each_img = to_tensor(each_img)[:3]
            current_img = torch.cat([current_img, each_img], dim=dim)
            diagonal = ((current_location[0]-start_location[0])**2 + (current_location[1]-start_location[1])**2)**(1/2)
            
        current_img = (to_img(current_img) if abs(ratio)<=1 else to_img(current_img.transpose(1,2)))
    elif distance_x < 0:
        ratio = distance_y/distance_x
        current_location = start_location
        # create biopy
        diagonal = 0
        current_img = (torch.zeros(3,pixel_radius,1) if abs(ratio)<=1 else torch.zeros(3,1,pixel_radius))
        size = ((1, pixel_radius) if abs(ratio)<=1 else (pixel_radius,1))
        dim = (2 if abs(ratio)<=1 else 1)
        move = (1 if y1 < y2 else -1)
        to_tensor = transforms.ToTensor()
        to_img = transforms.ToPILImage()
        while diagonal <= pixel_length:
            current_location = ((current_location[0] - 1, current_location[1] - ratio) if abs(ratio)<=1 else (current_location[0] - move/ratio, current_location[1] - move))
            each_img = wsi.read_region(location=(int(current_location[0]), int(current_location[1])),size=size, level=0)
            each_img = to_tensor(each_img)[:3]
            current_img = torch.cat([current_img, each_img], dim=dim)
            diagonal = ((current_location[0]-start_location[0])**2 + (current_location[1]-start_location[1])**2)**(1/2)
            
        current_img = (to_img(current_img) if abs(ratio)<=1 else to_img(current_img.transpose(1,2)))
    else:
        current_location = start_location
        diagonal = 0
        current_img = torch.zeros(3,1,pixel_radius)
        to_tensor = transforms.ToTensor()
        to_img = transforms.ToPILImage()
        move = (1 if y1 < y2 else -1)
        while diagonal <= pixel_length:
            current_location = (current_location[0], current_location[1] + move)
            each_img = wsi.read_region(location=(current_location[0], int(current_location[1])),size=(pixel_radius,1), level=0)
            each_img = to_tensor(each_img)[:3]
            current_img = torch.cat([current_img, each_img], dim=1)
            diagonal = ((current_location[0]-start_location[0])**2 + (current_location[1]-start_location[1])**2)**(1/2)
        current_img = to_img(current_img.transpose(1,2))
            
    return current_img

def make_biopsy(wsi_path, start_location, direct_location):
    wsi = openslide.OpenSlide(wsi_path)
    mpp_x = float(wsi.properties['openslide.mpp-x'])
    mpp_y = float(wsi.properties['openslide.mpp-y'])
    biopsy_length = np.random.normal(options['length'], options['length_sd']) * 1000 # convert to micron
    biopsy_radius = np.random.normal(options['radius'], options['radius_sd']) * 1000
    pixel_length = int(biopsy_length/mpp_x)
    pixel_radius = int(biopsy_radius/mpp_x)
    sample = get_biopsy(start_location=start_location, 
                 direct_location=direct_location, 
                 wsi=wsi, 
                 pixel_length=pixel_length, 
                 pixel_radius=pixel_radius)
    sample = sample.resize((int(pixel_length),int(pixel_radius)))
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