from PIL import Image
from math import *
import numpy

# Returns a metric for dissimilarity of each element in a list from its k nearest neighbours: sqrt of sum of squares of distances
def get_k_dist(lst, i, k):
        total = 0
        for j in xrange(1, k+1):
                if i-j > 0:
                        total += (lst[i] - lst[i-j])**2
                if i+j < len(lst):
                        total += (lst[i] - lst[i+j])**2
        return sqrt(total)

# Converts RGB coordinates to the Lab colour space (with XYZ space as an intermediate)
def rgb_to_lab(r, g, b):
    r, g, b = r/255.0, g/255.0, b/255.0
    if r > 0.04045:
        r = ((r + 0.055) / 1.055 )**2.4
    else:
        r /= 12.92

    if g > 0.04045:
        g = ((g + 0.055) / 1.055 )**2.4
    else:
        g /= 12.92
        
    if b > 0.04045:
        b = ((b + 0.055) / 1.055 )**2.4
    else:
        b /= 12.92
        
        # convert RGB to XYZ first, then XYZ to Lab
    x = (r*41.24 + g*35.76 + b*18.05) / 95.047
    y = (r*21.26 + g*71.52 + b*07.22) / 100
    z = (r*1.93 + g*11.92 + b*95.05) / 108.883
    if x > 0.008856:
        x **= (1.0/3.0)
    else:
        x = (7.787*x) + (16.0/116.0)
        
    if y > 0.008856:
        y **= (1.0/3.0)
    else:
        y = (7.787*y) + (16.0/116.0)
        
    if z > 0.008856:
        z **= (1.0/3.0)
    else:
        z = (7.787*z) + (16.0/116.0)
    return (116.0 * y) - 16, 500.0*(x - y), 200.0*(y - z)

# Calculates color difference using the CMC Delta_E formula
# Reference: http://en.wikipedia.org/wiki/Color_difference#CMC_l:c_.281984.29
def deltae_cmc(lab1, lab2):
    pl, pc = 2, 1 #2:1 ratio
    L1, a1, b1 = lab1
    L2, a2, b2 = lab2
    
    delta_L, delta_a, delta_b = L1 - L2, a1 - a2, b1 - b2
    C_1 = sqrt(a1**2 + b1**2)
    C_2 = sqrt(a2**2 + b2**2)
    H_1 = degrees(atan2(b1, a1))
    
    if H_1 < 0:
        H_1 = H_1 + 360
    
    F = sqrt(C_1**4 / (C_1**4 + 1900.0))
    if 164 <= H_1 and H_1 <= 345:
        T = 0.56 + abs(0.2 * cos(radians(H_1 + 168)))
    else:
        T = 0.36 + abs(0.4 * cos(radians(H_1 + 35)))
        
    if L1 < 16:
        S_L = 0.511
    else:
        S_L = (0.040975*L1) / (1 + 0.01765*L1)
        
    S_C = ((0.0638*C_1) / (1 + 0.0131*C_1)) + 0.638
    S_H = S_C * (F*T + 1 - F)    
    delta_C = C_1 - C_2
    try:
        delta_H = sqrt(delta_a**2 + delta_b**2 - delta_C**2)
    except ValueError:
        delta_H = 0.0
    
    L_group = delta_L / (pl * S_L)
    C_group = delta_C / (pc * S_C)
    H_group = delta_H / S_H    
    return sqrt(L_group**2 + C_group**2 + H_group**2)
    
# Call this function with the shredded image filepath to unshred.
# Default parameters: k=10 neighbours for k-distance calculation, std>2 standard deviations to define outliers in distance data, output file = glued.png
# If number of shreds is not provided, then glue() will attempt to autodetect the shreds.
def glue(fin, fout="glued.png", shredcount=None, k=10, stds=2):
    img = Image.open(fin)
    data = img.getdata()
    width, height = img.size

        # Get the pixel data for x,y coordinates and return the Lab coordinates, disregarding alpha
    def getpixel(x, y):
        pixel = data[y * width + x]
        return rgb_to_lab(pixel[0], pixel[1], pixel[2])

        # Return the total color difference of corresponding pixels in two columns of the image
    def coldiff(x1, x2):
        return sum([deltae_cmc(getpixel(x1, y), getpixel(x2, y)) for y in xrange(0, height)])

        # BONUS: Autodetect number of shreds by calculating successive column diffs and identifying outliers as shred boundaries
    if shredcount:
        shredwidth = width/shredcount
    else:
        print "autodetecting shreds"
        diffs = [coldiff(i-1, i) for i in xrange(1, width)]
        dists = [get_k_dist(diffs, i, k) for i in xrange(0, len(diffs))]
        threshold = numpy.mean(dists) + 2*numpy.std(dists)
        outliers = [i for i in xrange(0, width-1) if dists[i] > threshold]
        # Account for noise by only taking the most common width between shred boundaries
        widths = {}
        for i in xrange(1, len(outliers)):
                widths[outliers[i]-outliers[i-1]] = widths.get(i, 0) + 1
        shredwidth = max(widths.items(), key=lambda(tup): tup[1])[0]
        shredcount = width/shredwidth
        print "     Shreds:", shredcount, "\n     Width:", shredwidth, "px"
        
    shredbounds = [(i*shredwidth, (i+1)*shredwidth - 1) for i in xrange(0, shredcount)]
    newshreds = [shredbounds.pop()]
    print "Shred 1 glued"
        # Diff the left and right columns of each remaining shred with the working image, and glue in the best fit on the left or right
    while len(shredbounds):
        leftdiff, leftshred = min([(coldiff(shred[1], newshreds[0][0]), shred) for shred in shredbounds], key = lambda(tup): tup[0])
        rightdiff, rightshred = min([(coldiff(newshreds[-1][1], shred[0]), shred) for shred in shredbounds], key = lambda(tup): tup[0])
        if leftdiff < rightdiff:
            newshreds.insert(0, leftshred)
            shredbounds.remove(leftshred)
        else:
            newshreds.append(rightshred)
            shredbounds.remove(rightshred)
        print "Shred", len(newshreds), "glued"

        # Save the unshredded image to an output file
    glued = Image.new("RGBA", img.size)
    for i, shred in enumerate(newshreds):
        region = img.crop((shred[0], 0, shred[1]+1, height))
        glued.paste(region, (i*shredwidth, 0))
    glued.save(fout)
    print "Your 'secret' stuff is now in", fout+". Muahahahaha"
