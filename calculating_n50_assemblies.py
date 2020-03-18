#I am now heading to exercise two 
#I have a pseudocode at the moment which I would ideally like to implement in python
"""
while reading zipped assembly file line by line:
	if line matches "^>":
		word after that (e.g. RKZR02000001.1) = contig_id
		length_contig_id = 0
		append contig_id to contig_list
	else:
	#or alternatively elif line matches "^[ATGC]":
		split line by characters (nucleotides)
		length_contig_id += length of number of characters

for each element in contig_list:
	print contig_id , length_contig_id
	append length_contig_id to length_contig_list



#N50 calculation

sum_contig_lengths = 0

for each element in length_contig_list:
	sum_contig_lengths += length_contig_id

L50_length = sum_contig_lengths/2
temp_sum = 0
reverse sort length_contig_list (decreasing order, longest first)

for n (0..end_of_length_contig_list):
	temp_sum += length_contig_list[n]
	if (temp_sum >= L50_length):
		print "N50 = length_contig_list[n]"
"""

