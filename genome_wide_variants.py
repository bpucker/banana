### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.15 ###

### based on https://doi.org/10.1371/journal.pone.0164321 ###

__usage__ = """
					python genome_wide_variants.py
					--vcf <FULL_PATH_TO_INPUT_VCF>
					--out <FULL_PATH_TO_OUTPUT_DIR>
					
					optional:
					--res <INT, RESOLUTION>[1000000]
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys, os

# --- end of imports --- #


def load_variants_from_vcf( vcf_file ):
	"""! @brief loads the variant informaiton from a SnpEff output VCF file """
	
	snps_per_chr = {}
	indels_per_chr = {}
	
	tri_counter = 0
	
	with open( vcf_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if not "," in parts[4]:	#only biallelic variants
					if len( parts[3] ) == len( parts[4] ) and len( parts[3] ) == 1:
						try:
							snps_per_chr[ parts[0] ].append( parts[1] )
						except KeyError:
							snps_per_chr.update( { parts[0]: [ parts[1] ] } )
						
					elif len( parts[3] ) != len( parts[4] ):
						try:
							indels_per_chr[ parts[0] ].append( parts[1] )
						except KeyError:
							indels_per_chr.update( { parts[0]: [ parts[1] ] } )
				else:	#count triallelic variants
					tri_counter += 1
						
			line = f.readline()
	print "number of triallelic variants: " + str( tri_counter )
	
	return snps_per_chr, indels_per_chr


def generate_binned_values( lower_lim, upper_lim, chr_length, snps_per_chr, indels_per_chr, resolution ):
	"""! @brief group variants into bins """
	
	snp_data = []
	indel_data = []
	while True:
		if upper_lim >= chr_length:
			break
		else:
			snp_tmp = []
			indel_tmp = []
			for SNP in snps_per_chr:
				if SNP <= upper_lim and SNP > lower_lim:
					snp_tmp.append( 'X' )
			for indel in indels_per_chr:
				if indel <= upper_lim and indel > lower_lim:
					indel_tmp.append( 'X' )
			snp_data.append( len( snp_tmp ) )
			indel_data.append( len( indel_tmp ) )
		upper_lim += resolution
		lower_lim += resolution
	return max( snp_data ), max( indel_data ), snp_data, indel_data


def construct_plot( snps_per_chr_in, indels_per_chr_in, result_file, result_table, resolution ):
	"""! @brief construct variant over Col-0 genome distribution plot """
	
	# --- conversion of data into lists of lists --- #
	chr_lengths = []
	chr_names = []
	snps_per_chr = []
	indels_per_chr = []
	for key in sorted( snps_per_chr_in.keys() ):
		snps = snps_per_chr_in[ key ]
		indels = indels_per_chr_in[ key ]
		snps_per_chr.append( snps )
		indels_per_chr.append( indels )
		chr_lengths.append( max( snps+indels ) )
		chr_names.append( key )
	
	max_x_value = max( [ x for chro in snps_per_chr for x in chro ] + [ x for chro in indels_per_chr for x in chro ] )
	
	# --- generation of figure --- #
	fig, ax = plt.subplots()
	
	snp_scale = 1
	indel_scale = 1
	snp_data = []
	indel_data = []
	for idx, chr_length in enumerate( chr_lengths ):
		max_snp, max_indel, snp_temp, indel_temp = generate_binned_values( 0, resolution+0, chr_length, snps_per_chr[ idx ], indels_per_chr[ idx ], resolution )
		snp_data.append( snp_temp )
		indel_data.append( indel_temp )
		snp_scale = max( [ snp_scale, max_snp ] )
		indel_scale = max( [ indel_scale, max_indel ] )
	
	snp_scale = float( snp_scale )
	indel_scale = float( indel_scale )
	y_max = len( chr_lengths )
	
	ax2 = ax.twinx()
	
	with open( result_table, "w" ) as out:
		for idx, chr_length in enumerate( chr_lengths ):
			y = y_max-( idx*1.2 )
			x = resolution / 1000000.0
			
			ax.text( ( chr_length/ 1000000.0 ), y+0.3, chr_names[ idx ], ha="right" )
			
			# --- plotting SNP and InDel distribution --- #
			for i, snps in enumerate( snp_data[ idx ] ):
				indels = indel_data[ idx ][ i ]
				
				ax.plot( [ x*i+0.5*x, x*i+0.5*x ], [ y, y+ ( snps / snp_scale ) ], "-", color="lime" )
				ax2.plot( [ x*i+0.5*x, x*i+0.5*x ], [ y, y+ ( indels / indel_scale ) ], "-", color="magenta" )
			
			ax.plot( [ 0, 0 ], [ y, y+1 ], color="black" )
			ax.text( 0, y+1, str( int( snp_scale ) ), ha="right", fontsize=5 )
			ax.text( 0, y+0.5, str( int( snp_scale / 2 ) ), ha="right", fontsize=5 )
			ax.text( 0, y, "0", ha="right", fontsize=5 )
			
			ax.plot( [ max_x_value, max_x_value ], [ y, y+1 ], color="black" )
			ax.text( max_x_value, y+1, str( int( indel_scale ) ), ha="right", fontsize=5 )
			ax.text( max_x_value, y+0.5, str( int( indel_scale / 2 ) ), ha="right", fontsize=5 )
			ax.text( max_x_value, y, "0", ha="right", fontsize=5 )
			
			# --- writing data into output table --- #
			out.write( 'Chr' + str( idx+1 ) + "SNVs:\t" + '\t'.join( map( str, snp_data ) ) + '\n' )
			out.write( 'Chr' + str( idx+1 ) + "InDels:\t" + '\t'.join( map( str, indel_data ) ) + '\n' )
	
	ax.set_xlabel( "genomic position [ Mbp ]" )
	ax.set_ylabel( "number of SNVs per interval" )
	ax2.set_ylabel( "number of InDels per interval" )
	
	ax.set_xlim( 0, max( chr_lengths ) / 1000000.0 )	
	
	ax.legend( handles=[ mpatches.Patch(color='magenta', label='InDels'), mpatches.Patch(color='lime', label='SNVs') ], prop={'size':10} )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	
	ax.get_yaxis().set_ticks([])
	ax2.get_yaxis().set_ticks([])
	ax.yaxis.labelpad = 15
	ax2.yaxis.labelpad = 15
	
	
	plt.subplots_adjust( left=0.1, right=0.9, top=0.99, bottom=0.1 )
	fig.savefig( result_file, dpi=600 )
	
	plt.close('all')


def main( arguments ):
	"""! @brief run everyting """
	
	vcf_file = arguments[ arguments.index( '--vcf' )+1 ]
	output_dir = arguments[ arguments.index( '--out' )+1 ]
	
	if '--res' in arguments:
		resolution = int( arguments[ arguments.index( '--res' )+1 ] )
	else:
		resolution = 1000000
	
	if output_dir[ -1 ] != "/":
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	result_file = output_dir + "genome_wide_small_variants.png"
	result_table = output_dir + "genome_wide_small_variants.txt"
	
	snps_per_chr, indels_per_chr = load_variants_from_vcf( vcf_file )
	
	print "number of SNVs: " + str( len( [ x for each in snps_per_chr.values() for x in each ] ) )
	print "number of InDels: " + str( len( [ x for each in indels_per_chr.values() for x in each ]) )
	
	construct_plot( snps_per_chr, indels_per_chr, result_file, result_table, resolution )


if '--vcf' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
