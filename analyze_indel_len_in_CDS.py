### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python analyze_indel_len_CDS.py
					--vcf <FULL_PATH_TO_VCF_FILE>
					--gff <FULL_PATH_TO_GFF3_FILE>
					--out <FULL_PATH_TO_OUTPUT_FOLDER>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""


import os, sys, math
import matplotlib.pyplot as plt

# --- end of imports --- #

def load_all_CDS_positions( gff ):
	"""! @brief load all CDS positions """
	
	CDS_pos = {}
	with open( gff, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "CDS":
					start, end = map( int, parts[3:5] )
					while start < end:
						key = parts[0] + "_%_" + ( str( start ).zfill( 8 ) )
						try:
							CDS_pos[ key ]
						except KeyError:
							CDS_pos.update( { key: None } )
						start += 1
			line = f.readline()
	return CDS_pos


def generate_figure( indel_lengths, figfile ):
	"""! @brief generate boxplot with InDel lengths """
	
	# --- check InDel lengths --- #
	len_to_check = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 ]
	values = []
	for length in len_to_check:
		try:
			values.append( math.log10( indel_lengths.count( length ) ) )
		except ValueError:
			values.append( 0 )
			
	
	# --- generate figure --- #	
	fig, ax = plt.subplots()
	ax.bar( len_to_check, values, color="lime" )
	ax.set_xlabel( "InDel length distribution" )
	ax.set_ylabel( "log10( number of InDels )" )
	
	ax.set_xlim( 0, 30.5 )
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	plt.subplots_adjust( left=0.1, right=0.99, top=0.99, bottom=0.1 )
	
	fig.savefig( figfile, dpi=300 )
	plt.close( "all" )


def main( arguments ):
	"""! @brief run everything """
	
	vcf = arguments[ arguments.index('--vcf')+1 ]
	gff = arguments[ arguments.index('--gff')+1 ]
	output_dir = arguments[ arguments.index('--out')+1 ]
	
	if not output_dir[-1] == "/":
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	CDS_pos = load_all_CDS_positions( gff )
	
	
	CDS_indel_lens = []
	other_indel_lens = []
	
	with open( vcf, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if not "," in parts[1]:
					if len( parts[3] ) != len( parts[4] ):
						try:
							CDS_pos[ parts[0] + "_%_" + ( parts[1].zfill( 8 ) ) ]
							CDS_indel_lens.append( abs( len( parts[3] ) - len( parts[4] ) ) )
						except KeyError:
							other_indel_lens.append( abs( len( parts[3] ) - len( parts[4] ) ) )
			line = f.readline()
	
	print "number of InDels in CDS: " + str( len( CDS_indel_lens ) )
	print "number of InDels outside CDS: " + str( len( other_indel_lens ) )
	print "total CDS length: " + str( len( CDS_pos.keys() ) )
	
	CDS_file = output_dir + "CDS_InDel_lengths.png"
	other_file = output_dir + "other_InDel_lengths.png"
	
	generate_figure( CDS_indel_lens, CDS_file )
	generate_figure( other_indel_lens, other_file )


if '--vcf' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
