### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python ploidy_check.py
					--vcf <FULL_PATH_TO_INPUT_VCF>
					--out <FULL_PATH_TO_OUTPUT_DIRECTORY>
					"""

import matplotlib.pyplot as plt
import sys, os


# --- end of imports --- #

def main( arguments ):
	"""! @brief run everything """
	
	vcf_file = arguments[ arguments.index( '--vcf' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	if output_dir[ -1 ] != "/":
		output_dir += "/"
	
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	fig_file = output_dir + "allele_frequencies.png"
	cov_fig_file = output_dir + "variant_coverages.png"

	values = []
	coverage = []


	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[-1][:3] == "0/1":
					x, y = map( float, parts[-1].split(':')[1].split(',')[:2] )
					coverage.append( x+y )
					if x+y > 20:
						values.append( x / ( x+y ) )
				elif parts[-1][:3] == "1/1":
					x, y = map( float, parts[-1].split(':')[1].split(',')[:2] )
					coverage.append( x+y )
			line = f.readline()

	# --- generation of variant coverage histogram --- #
	fig, ax = plt.subplots()
	ax.hist( coverage, bins=10000, color="lime" )
	ax.set_xlim( 0, 400 )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	plt.subplots_adjust( left=0.13, right=0.97, top=0.99, bottom=0.1 )
	ax.set_xlabel( "sequencing coverage" )
	ax.set_ylabel( "number of variants" )
	
	fig.savefig( cov_fig_file, dpi=300 )

	
	# --- generation of variant frequency histogram --- #
	fig, ax = plt.subplots()
	ax.hist( values, bins=100, color="lime" )
	ax.set_xlabel( "allele frequency" )
	ax.set_ylabel( "number of variants" )
	
	ax.plot( [ 0.33, 0.33 ], [ 0, 100000 ], color="black" )
	ax.plot( [ 0.5, 0.5 ], [ 0, 100000 ], color="black" )
	ax.plot( [ 0.66, 0.66 ], [ 0, 100000 ], color="black" )
	
	ax.set_xlim( 0, 1 )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	plt.subplots_adjust( left=0.15, right=0.99, top=0.99, bottom=0.1 )
	
	fig.savefig( fig_file, dpi=300 )
	
	plt.close( "all" )


if '--vcf' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
