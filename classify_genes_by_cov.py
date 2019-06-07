### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.15 ###

__usage__ = """
					python classify_genes_by_cov.py
					--gff <FULL_PATH_TO_GFF3_FILE>
					--cov <FULL_PATH_TO_COVERAGE_FILE>
					--out <FULL_PATH_TO_OUTPUT_DIR>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import matplotlib.pyplot as plt
import numpy as np
import sys, os


def load_cov( cov_file ):
	"""! @brief load all information from coverage file """
	
	cov = {}
	with open( cov_file, "r" ) as f:
		line = f.readline()
		header = line.split('\t')[0]
		tmp = []
		while line:
			parts = line.strip().split('\t')
			if parts[0] != header:
				cov.update( { header: tmp } )
				header = parts[0]
				tmp = []
			tmp.append( float( parts[-1] ) )
			line = f.readline()
		cov.update( { header: tmp } )
	return cov


def load_gene_positions( gff3_file ):
	"""! @brief load all gene positions from given gff3 file """
	
	gene_pos = {}
	with open( gff3_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				if parts[2] == "gene":
					# --- get gene ID --- #
					if ";" in parts[-1]:
						ID = parts[-1].split(';')[0].split('ID=')[-1]
					else:
						ID = parts[-1].split('ID=')[-1]
					gene_pos.update( { ID: { 'chr': parts[0], 'start': int( parts[3] ), 'end': int( parts[4] ) } } )
			line = f.readline()
	return gene_pos


def get_gene_covs( cov, gene_pos ):
	"""! @brief get average coverage per gene """
	
	gene_covs = {}
	for gene in gene_pos.keys():
		try:
			gene_covs.update( { gene: np.mean( cov[ gene_pos[ gene ]['chr'] ][ gene_pos[ gene ]['start']:gene_pos[ gene ]['end'] ] ) } )
		except KeyError:
			print gene
	return gene_covs


def count_chromosomes( gene_pos ):
	"""! @brief count number of chromosomes to get proper height """
	
	chromosomes = []
	for gene in gene_pos:
		chromosomes.append( gene_pos[ gene ]['chr'] )
	chr_names = sorted( list( set( chromosomes ) ) )
	return len( chr_names ) + 1, chr_names


def generate_heatmap( gene_pos, gene_covs, fig_file ):
	"""! @brief construct heatmap """
	
	upper_cutoff = 300.0
	lower_cutoff = 0.0
	
	fig, ax = plt.subplots( figsize=( 8, 4 ) )
	
	y_max, chr_names = count_chromosomes( gene_pos )	#number of chromosomes + 1
	x_values = []
	y_values = []
	color_values = []
	
	for gene in gene_pos:
		x_values.append( gene_pos[ gene ]['start'] / 1000000.0 )
		y_values.append( y_max - chr_names.index( gene_pos[ gene ]['chr'] ) - 1 )
		c = max( [ min( [ gene_covs[ gene ], upper_cutoff ] ), lower_cutoff ] )
		color_values.append( c )
	
	ax.scatter( x_values, y_values, c=color_values, cmap="cool", s=1, edgecolors="none" )
	
	ax.set_xlim( 0, max( x_values ) )
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.get_yaxis().set_ticks([])
	
	legend_x = [ max( x_values ) - ( 0.1 * max( x_values ) ), max( x_values ) - ( 0.1 * max( x_values ) ), max( x_values ) - ( 0.1 * max( x_values ) ), max( x_values ) - ( 0.1 * max( x_values ) ) ]
	legend_y = [ y_max-2, y_max-3, y_max-4, y_max-5 ]
	legend_color = [ 0.0, 100.0, 200.0, 300.0  ]
	ax.scatter( legend_x, legend_y, c=legend_color, cmap="cool", s=10, edgecolors="none" )
	
	ax.text( max( x_values ) - ( 0.08 * max( x_values ) ), y_max-2, "0x", fontsize=5, ha="left", va="center" )
	ax.text( max( x_values ) - ( 0.08 * max( x_values ) ), y_max-3, "100x", fontsize=5, ha="left", va="center" )
	ax.text( max( x_values ) - ( 0.08 * max( x_values ) ), y_max-4, "200x", fontsize=5, ha="left", va="center" )
	ax.text( max( x_values ) - ( 0.08 * max( x_values ) ), y_max-5, "300x", fontsize=5, ha="left", va="center" )
	
	ax.set_xlabel( "position on chromosome [ Mbp ]" )
	
	plt.subplots_adjust( left=0.01, right=0.99, top=0.99, bottom=0.25 )
	
	fig.savefig( fig_file, dpi=300 )


def main( arguments ):
	"""! @brief run everything """
	
	gff3_file = arguments[ arguments.index( '--gff' ) +1 ]
	cov_file = arguments[ arguments.index( '--cov' ) +1 ]
	output_dir = arguments[ arguments.index( '--out' ) +1 ]
	
	if not output_dir[ -1 ] == "/":
		output_dir += "/"
	if not os.path.exists( output_dir ):
		os.makedirs( output_dir )
	
	fig_file = output_dir + "gene_coverage_heatmap.png"
	
	gene_pos = load_gene_positions( gff3_file )
	
	cov = load_cov( cov_file )
	
	gene_cov = get_gene_covs( cov, gene_pos )
	
	generate_heatmap( gene_pos, gene_cov, fig_file )


if '--gff' in sys.argv and '--cov' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
