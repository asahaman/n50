import glob
import re
import gzip
import operator

file1 = open("summary_concise.txt","w")
file2 = open("summary_verbose.txt","w")

#intializing assembly count
assembly_count = 0

#iterating through each assembly file
for file in glob.glob('*.gz'):
    file2.write ("Opened file {}\n".format(file))
    assembly_count += 1
    with gzip.open(file, mode = 'rt') as f:

        #initializing contig length dictionary and counting contigs
        contig_length_dict = {}
        contig_count = 0

        #reading assembly file line by line
        for line in f:

            #capturing assembly name. Although this gets updated with every contig header,
            #it remains same for every file due to the nature of the regex pattern
            assembly = re.findall(r'^>([A-Z]+?0\d).+',line)
            if len(assembly) > 0:
                assembly_name = assembly[0] 

            #capturing contig name, initializing contig length to zero, incrementing contig count       
            x = re.findall(r'^>(.+?)\s', line)
            if len(x) > 0:
              contig_name = x[0]
              contig_length = 0
              contig_count += 1
            
            #if the program doesn't find a contig header, the seqeunce length of every line gets added
            #and the contig length gets updated within the contig length dictionary
            else:
              contig_length += len(line) - 1
              contig_length_dict[contig_name] = contig_length

        file2.write("File has assembly named {} with following contigs and respective lengths:\n".format(assembly_name))
        
        #The assembly length is calculated by summing all contig lengths
        #The length of every contig from the contig length dictionary is written to a file
        total_assembly_len = 0
        for key in contig_length_dict:
            file2.write("\t{}\t{}\n".format(key,contig_length_dict[key]))
            total_assembly_len += contig_length_dict[key]
        file2.write("Assembly {} has {} contigs\n\n".format(assembly_name,contig_count))
        
        #L50 length is calculated
        l50_len = total_assembly_len/2
        
        #contig length dictionary is sorted by values to an ordered list of tuples
        rev_sorted_contig_length_dict = sorted(contig_length_dict.items(), key=operator.itemgetter(1),reverse=True)
        temp_sum = 0
        
        #iterate through decreasing contig lengths and breaking the loop once N50 length is calculated
        for item in rev_sorted_contig_length_dict:
            temp_sum += item[1]
            if temp_sum >= l50_len: 
                file1.write("N50 length for assembly {} is {}\n".format(assembly_name, item[1]))
                break