#Inspired from Liang Tong's bargraph. Pygraph written in python 2.7.3 puts protein sequence information
#(alignment, secondary structure prediction, disorder prediction, domain organization)
#into one graphics file.

#Instructions: 
#Edit lines following "Enter".
#To excecute, type at the command line "python pygraph_final.py"
#postscript file out.ps will write in your working directory

#Make sure to have the graphics library file 'http://mcsp.wartburg.edu/zelle/python/graphics.py'
#and relevent text files generated from clustalw2 and UCL 'http://bioinf.cs.ucl.ac.uk'
#in your working directory.
#If any file is not provided, comment-out the draw function at the bottom of the page.

from graphics import *

WIDTH = 1000
HEIGHT = 250        
win = GraphWin("pygraph", WIDTH, HEIGHT)

#Enter number of sequences in clustal alignment file: i.e 7 if seven sequences used for aligning. 
seqs = 7

#Enter the target sequence from the top of alignment file: i.e 5 if 5th from top 
target = 1

#Enter the clustalw output in text format: i.e aug8_list_clustalw.
infile = open("cin1_clustalw" , "r")
line = infile.readlines()
infile.close()

#Enter boundry of sequence alignment: ie. 17th - 76th column in the clustal text. (varies  due to name length)
bound = [17, 76]

#Enter the psipred output in text format: i.e sc_cin4_AA
infile = open("sc_cin1_AA" , "r")
line2 = infile.readlines()
infile.close()

#Enter the disopred output in text format: i.e sc_cin4_DISOPRED
infile = open("sc-cin1_DISOPRED" , "r")
line3 = infile.readlines()
infile.close()

#Enter a name for target : i.e sc_cin4
NAME = 'sc_cin1'

#Enter domain range in paranthesis separated by commas. ie. (1,10), (20,45) 
domains = ()

#read clustalw alignment
def extract_clustal(seqs, target, bound, line):

    n = seqs + 2

    # extracts conservation symbols( *, :, . ) from clustal alignment file.     
    clustal = []
    for i in range (seqs+3,len(line),n):   #extracts "*:.." symbols from clustal file into a list
        clustal.append(line[i])

    clustal1 = []        
    for i in clustal:           #removes the extra space in front of each line
        a = i[bound[0]-1:bound[1]]
        clustal1.append(a)
    

    clustalw = []
    for i in clustal1:          #puts conservation symbols into a list "clustalw". one for each AA
        clustalw += i 
    
    #remove last '\n'
    del clustalw[-1]
        
    # removes gaps between symbols due to gaps in the sequence.  
    seq0 = []
    for i in range (target+2, len(line), n):      #extract target sequence from clustal file into a list
        seq0.append(line[i])
                
    seq1 = []
    for i in seq0:                          #removes the extra space in fron of each line
        a = i[bound[0]-1:bound[1]]
        seq1.append(a)

    seq = []
    for i in seq1:
        seq += i                       #puts target sequence from algn file into a list "seq"


    #remove \n and integers from seq list. 
    while len(seq) != len(clustalw):
        del seq[-1]
    
    for i in range (len(seq)):
        if seq[i] == '-' and clustalw[i] == ' ':    #tags the gaps from target algnined sequence
            seq[i] = '$'
            clustalw[i] = '$'

    
    while '$' in seq:                           #removes gaps from sequence
        seq.remove('$')
    while '$' in clustalw:                      #removes gaps from the aligned symbols "*:." 
        clustalw.remove('$')

    return clustalw
    
    

#read secondary structure probality and symbols from PSIPRED text file. 
def extract_secondary(line):

    # extracts probability of secondary structures and put into list.
    psi = []
    for i in range (2,len(line),6):
        psi.append(line[i])            #puts 2ndary probability from each line in a list
    psi1 = []        
    for i in psi:                       #removes the extra space in front of each line
        a = i[6:]
        psi1.append(a)
    psi2 = []
    for i in psi1:                      #removes the "\n" at the end of each line
        b = i.rstrip('\n')
        psi2.append(b)
    psipred = []
    for i in psi2:                          #concatenate probabilities from list in one character
        psipred += i + ''
    psipred = [int(i) for i in psipred]     #convert number strings to number integers


    # extracts symbol for secondary structure C= coil, E=strand, H=helix    
    C = []
    for i in range (3, len(line), 6):        #puts 'CEH' symbols from each line2 in a list
        C.append(line[i])
    E = []
    for i in C:                             #removes the extra space in front of each line
        a = i[6:]
        E.append(a)
    H = []              
    for i in E:                             #removes the "\n" at the end of each line
        b = i.rstrip('\n')
        H.append(b)
    CEH = []
    for i in H:                             #concatenate probabilities from list in one character
        CEH += i + ''


    return psipred, CEH
    

#read DISOPRED file
def extract_disorder(line):

    # extract disorder probability and put into list.
    dis = []
    for i in range (2,len(line) -4 ,5):   #puts disorder probability from each line in a list
        dis.append(line[i])            
    diso = []
    for i in dis:                         #removes the extra space in front of each line
        a = i[6:]
        diso.append(a)
    disop = []              
    for i in diso:                         #removes the "\n" at the end of each line
        b = i.rstrip('\n')
        disop.append(b)
    disopred = []
    for i in disop:                         #concatenate probabilities from list in a list
        disopred += i + ''
    disopred = [int(i) for i in disopred]

    return disopred

           

def draw_clustal(n):


    # color-code low-to-high conservation: black-blue-yellow-red
    p1 = Point ((0.005*WIDTH*0)+0.01*WIDTH, 85)
    p2 = Point ((0.005*WIDTH*0)+0.01*WIDTH, 105)
    rec = Line (p1, p2)
    rec.setFill('black')
    rec.setWidth(0.005*WIDTH)
    rec.draw(win)

    p1 = Point ((0.005*WIDTH*1)+0.01*WIDTH, 85)
    p2 = Point ((0.005*WIDTH*1)+0.01*WIDTH, 105)
    rec = Line (p1, p2)
    rec.setFill('blue')
    rec.setWidth(0.005*WIDTH)
    rec.draw(win)

    p1 = Point ((0.005*WIDTH*2)+0.01*WIDTH, 85)
    p2 = Point ((0.005*WIDTH*2)+0.01*WIDTH, 105)
    rec = Line (p1, p2)
    rec.setFill('yellow')
    rec.setWidth(0.005*WIDTH)
    rec.draw(win)

    p1 = Point ((0.005*WIDTH*3)+0.01*WIDTH, 85)
    p2 = Point ((0.005*WIDTH*3)+0.01*WIDTH, 105)
    rec = Line (p1, p2)
    rec.setFill('red')
    rec.setWidth(0.005*WIDTH)
    rec.draw(win)

    
    
    # clustal sequence alignment color bars
    count = 1
    for i in n:
        p1 = Point ((((0.95*WIDTH)/len(n))*count)+0.05*WIDTH, 85)
        p2 = Point ((((0.95*WIDTH)/len(n))*count)+0.05*WIDTH, 105)

        if i == '*':
            rec = Line (p1, p2)
            rec.setFill('red')
        elif i == ':':
            rec = Line (p1, p2)
            rec.setFill('yellow')
        elif i == '.':
            rec = Line (p1, p2)
            rec.setFill('blue')
        else:
            rec = Line (p1, p2)
            rec.setFill('black')
        rec.setWidth(((0.94*WIDTH)/len(n)))
        rec.draw(win)
        count += 1

    # tik-mark every 50 residues
    count = 1
    for i in range(len(n)):
        if count%50 == 0:
            p1 = Point ((((0.95*WIDTH)/len(n))*count)+0.05*WIDTH, 75)
            p2 = Point ((((0.95*WIDTH)/len(n))*count)+0.05*WIDTH, 80)
            rec = Line (p1,p2)
            rec.setFill('black')
            rec.setWidth(((0.94*WIDTH)/len(n)))
            rec.draw(win)
        count += 1

    #draw tik-mark      
    text = Text(Point((((0.95*WIDTH)/len(n))*50)+0.05*WIDTH, 65), '50')
    text.setSize(10)
    text.setFace('helvetica')
    text.draw(win)
    
    

def draw_secondary(sec_prob, sec_sym):

    # color code seconday structure. (to use threshold uncomment/add to condition)
    count = 1
    for i in range(len(sec_sym)):
        p1 = Point ((((0.95*WIDTH)/len(sec_sym))*count)+0.05*WIDTH, 110)
        p2 = Point ((((0.95*WIDTH)/len(sec_sym))*count)+0.05*WIDTH, 130)

        if sec_sym[count-1] == 'H':       #sec_prob[count] > 5 and 
            rec = Line(p1,p2)
            rec.setFill('cyan')

        elif sec_sym[count-1] == 'E':     #sec_prob[count] > 5 and 
            rec = Line(p1, p2)
            rec.setFill('magenta')

        else:
            rec = Line(p1,p2)
            rec.setFill('gray')

        rec.setWidth(((0.94*WIDTH)/len(sec_sym)))
        rec.draw(win)
        count += 1

    #draw legend text
    text = Text(Point(0.02*WIDTH, 115), 'alpha')
    text.setSize(10)
    text.setFill('cyan')
    text.setFace('helvetica')
    text.draw(win)

    text = Text(Point(0.02*WIDTH, 125), 'beta')
    text.setSize(10)
    text.setFill('magenta')
    text.setFace('helvetica')
    text.draw(win)
        
        
def draw_disorder(diso_prob):
            
    #color code disorder prediction. = or > 8 is disordered
    count = 1
    for i in range(len(diso_prob)):
        p1 = Point ((((0.95*WIDTH)/len(diso_prob))*count)+0.05*WIDTH, 135)
        p2 = Point ((((0.95*WIDTH)/len(diso_prob))*count)+0.05*WIDTH, 155)

        if diso_prob[count-1] >= 8:
            rec = Line (p1,p2)
            rec.setFill('black')

        else:
            rec = Line(p1,p2)
            rec.setFill('gray')

        rec.setWidth(((0.94*WIDTH)/len(diso_prob)))
        rec.draw(win)
        count += 1

    #draw legend text       
    text = Text(Point(0.025*WIDTH, 145), 'disorder')
    text.setSize(10)
    text.setFill('black')
    text.setFace('helvetica')
    text.draw(win)


def draw_domain(domains, n):

    #put domain borders in ranges
    ranges = []
    for i in domains:
        a = range(i[0], i[1]+1)
        ranges.append(a)

    #concatenate ranges in a list  
    count = []
    for i in ranges:
        count += i

    #color domain ranges.     
    for i in count:
        p1 = Point ((((0.95*WIDTH)/len(n))*i)+0.05*WIDTH, 160)
        p2 = Point ((((0.95*WIDTH)/len(n))*i)+0.05*WIDTH, 180)
        rec = Line(p1,p2)
        rec.setFill(color_rgb(0,100,0))
        rec.setWidth(((0.94*WIDTH)/len(n)))
        rec.draw(win)

    #put entire length of protein in a range
    ranges_n = range(1, len(n)+1)
    
    #subtract domains 
    for i in count:
        if i in ranges_n:
            ranges_n.remove(i)

    #color the rest of the protein gray
    for i in ranges_n:
        p1 = Point ((((0.95*WIDTH)/len(n))*i)+0.05*WIDTH, 160)
        p2 = Point ((((0.95*WIDTH)/len(n))*i)+0.05*WIDTH, 180)
        rec = Line(p1,p2)
        rec.setFill('gray')
        rec.setWidth(((0.94*WIDTH)/len(n)))
        rec.draw(win)
        
        
    
    #draw legend text       
    text = Text(Point(0.025*WIDTH, 170), 'domains')
    text.setSize(10)
    text.setFill(color_rgb(0,100,0))
    text.setFace('helvetica')
    text.draw(win)
    
    


draw_clustal(extract_clustal(seqs, target, bound, line))
draw_secondary(extract_secondary(line2)[0], extract_secondary(line2)[1])
draw_disorder(extract_disorder(line3))
draw_domain(domains, extract_clustal(seqs, target, bound, line))

text = Text(Point(0.075*WIDTH,45), NAME)
text.setSize(10)
text.setFace('helvetica')
text.draw(win)


win.postscript(file="out.ps" , colormode='color')

text = Text(Point(0.45*WIDTH, 0.95*HEIGHT), 'out.ps written, click screen to quit')
text.setSize(20)
text.setFace('helvetica')
text.draw(win)

win.getMouse()
win.close()

   
