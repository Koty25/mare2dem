# This parser works for the following Mare2dem files:
# .resistivity
# .topo
# .poly and its several parts
# More about each file and its segments can be found at:
# https://mare2dem.bitbucket.io/file_formats.html
#
# In order to properly generate the file which contains the
# resistivity and the regions together, the method parse_resistivity
# calls the method parse_poly and generates the needed "_regional_attr.csv" file.
# With both the resistivity and the regional_attr parsed, then we can merge 
# the two information.
#
# One can still parse only the poly file by simply executing the parser
# providing the desired .poly as input.
#
# Run the Script with:
# python3 parse_input.py -input <filename>

import argparse
import os
import re 

def parse_topo(infile):
    outfile = infile + ".csv"
    with open(infile) as fin, open(outfile, 'w') as fout:
        for line in fin:
            fout.write(line.replace('\t', ','))
    
    fin.close()
    fout.close()


def parse_poly(infile):
    count = 0
    with open(infile) as fin:
        while count < 4:
            if count is 0:
                header = ["vertex_id","y","z","attribute","boundary_marker"] 
                outfile = infile + "_vertices.csv"
            elif count is 1:
                header = ["edge_id", "start_point", "end_point", "boundary_marker"]
                outfile = infile + "_segments.csv"
            elif count is 2: 
                header = ["hole_id","y","z"] 
                outfile = infile + "_holes.csv"
            else:
                header = ["region_id","y","z","attribute","maximum_area"] 
                outfile = infile + "_regional_attr.csv"

            num = fin.readline().split()[0]
            if int(num) != 0:
                with open(outfile, 'w') as fout:
                    fout.write(",".join(header)+"\n")
                    for i in range(int(num)):
                        #splitting to remove the whitespace at the end of the line
                        values=fin.readline().split(" ")                        
                        fout.write(",".join(values[:len(values)-1])+"\n")

            count += 1
    
    fin.close()
    fout.close()


def parse_resistivity(infile):
    # Here we just parse the resistivity file, but we still need to 
    # parse the poly information:
    polyf = re.sub('\d+.resistivity', 'poly', infile)
    parse_poly(polyf)
    # now we want the ".poly_regional_attr.csv" file
    polyf = re.sub('poly', 'poly_regional_attr.csv', polyf)
    
    outfile = infile + ".csv"
    with open(infile) as fin, open(polyf) as fpoly, open(outfile, 'w') as fout:
        header = "index,rho,y,z\n"
        fout.write(header)

        line = next(fin)
        while not line.startswith((' Number of regions:','Number of regions:')):
            line = next(fin)

        # just skipping the header
        line = next(fin)
        poly = next(fpoly)

        # some resistivity files have a different amount of spaces
        # for which the previous Script didn't work
        # so this regex fixes it and works for every type
        regex = '\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)'
        list_data = re.findall(regex, fin.read())

        for data in list_data:
            poly = fpoly.readline()
            new_line = ','.join(data[:2])
            # this weird line selects the y and z values from the attributes file,
            # ignoring the index because we already have it
            poly = ','.join(poly.split(',')[1:3])
            new = new_line + ',' + poly

            fout.write(new)
            fout.write('\n')
    
    fin.close()
    fpoly.close()
    fout.close()

def main():
    parser = argparse.ArgumentParser(description="Script that parses Mare2dem input files into csv")
    parser.add_argument('-input', help='Location of the input file', required=True)
    
    args = parser.parse_args()
    i_file = args.input

    if i_file is None:
        print("Inform the input file:\npython3 parse_input.py -input <file_name>")
        exit()

    filename, extension = os.path.splitext(i_file)

    if extension == ".topo":
        parse_topo(i_file)
    elif extension == ".poly":
        parse_poly(i_file)
    elif extension == ".resistivity":
        parse_resistivity(i_file)
    else:
        print("Unkown input extension")

main()
