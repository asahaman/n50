import glob
import re
import gzip

assembly_count = 0
for file in glob.glob('*.gz'):
    print ("Opened file {}".format(file))
    assembly_count += 1
    with gzip.open(file, mode = 'rt') as f:
        contig_length_dict = {}
        contig_count = 0
        for line in f:
            assembly = re.findall(r'^>([A-Z]+?0\d).+',line)
            if len(assembly) > 0:
                assembly_name = assembly[0]    
            x = re.findall(r'^>(.+?)\s', line)
            if len(x) > 0:
              contig_name = x[0]
              contig_length = 0
              contig_count += 1
            else:
              contig_length += len(line) - 1
              contig_length_dict[contig_name] = contig_length
        print("File has assembly named {} with following contigs and respective lengths:".format(assembly_name))
        for key in contig_length_dict:
            print("\t",key,"\t",contig_length_dict[key])
        print("Assembly {} has {} contigs\n".format(assembly_name,contig_count))    

