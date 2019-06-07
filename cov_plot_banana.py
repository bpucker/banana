### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.1 ###

__usage__ = """
					python cov_plot_banana.py
					--in <FULL_PATH_TO_COVERAGE_FILE>
					--out <FULL_PATH_TO_OUTPUT_FILE>
					
					--res <RESOLUTION, WINDOW_SIZE_FOR_COVERAGE_CALCULATION>
					--sat <SATURATION, CUTOFF_FOR_MAX_COVERAGE_VALUE>
					"""

import sys, os
import matplotlib.pyplot as plt
import numpy as np

# --- end of imports --- #


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


def generate_plot( cov, out_file, resolution, saturation ):
	"""! @brief generate figure """
	
	fig, ax = plt.subplots( figsize=( 10, 7 ) )
	
	ymax = 12	#len( cov.keys() )+1
	max_value = 0
	collected_values = {}
	
	# --- generate list for plotting --- #
	for idx, key in enumerate( sorted( cov.keys() )[:12] ):
		y = ymax-idx-1
		x = []
		blocks = [ cov[ key ] [ i : i + resolution ] for i in xrange( 0, len( cov[ key ] ), resolution ) ]
		for block in blocks:
			x.append( min( [ np.mean( block ), saturation ] ) )
		max_value = max( [ max_value, max( x ) ] )
		collected_values.update( { key: x } )
	
	# --- plot values --- #
	max_value = float( min( [ saturation, max_value ] ) )
	for idx, key in enumerate( sorted( cov.keys() )[:12] ):
		y = ymax - ( idx*1.3 )
		x = []
		for each in collected_values[ key ]:
			x.append( y + min( [ 1, ( each / max_value ) ] ) )
		
		ax.scatter( np.arange( 0, len( x ), 1 ), x, s=1, color="lime" )
		
		ax.plot( [ 0, len( x ) ], [ y+( 0 / max_value ), y+( 0 / max_value ) ], color="black" , linewidth=0.1)
		ax.plot( [ 0, len( x ) ], [ y+( 50 / max_value ), y+( 50 / max_value ) ], color="black" , linewidth=0.1)
		ax.plot( [ 0, len( x ) ], [ y+( 100 / max_value ), y+( 100 / max_value ) ], color="black" , linewidth=0.1)
		ax.plot( [ 0, len( x ) ], [ y+( 150 / max_value ), y+( 150 / max_value ) ], color="black" , linewidth=0.1)
		ax.plot( [ 0, len( x ) ], [ y+( 200 / max_value ), y+( 200 / max_value ) ], color="black" , linewidth=0.1)
		
		ax.plot( [ 0, 0 ], [ y, y+1 ], color="black", linewidth=1, markersize=1 )
		ax.text( 0, y+1, str( int( max_value ) ), ha="right", fontsize=5 )
		ax.text( 0, y+0.5, str( int( max_value / 2 ) ), ha="right", fontsize=5 )
		ax.text( 0, y, "0", ha="right", fontsize=5 )
		
	ax.set_xlabel( "position on chromosome [ Mbp ]" )
	ax.set_ylabel( "coverage" )
	
	ax.set_xlim( 0, 4800 )
	
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['left'].set_visible(False)
	ax.get_yaxis().set_ticks([])
	ax.yaxis.labelpad = 10
	
	ax.xaxis.set_ticks( np.arange( 0, 4900, 100 ) )
	labels = map( str, np.arange( 0, 49, 1 ) )
	ax.set_xticklabels( labels )
	
	plt.subplots_adjust( left=0.03, right=0.999, top=0.99, bottom=0.1 )
	
	fig.savefig( out_file, dpi=300 )
	plt.close( "all" )


def generate_hist( cov_values, outputfile, saturation, resolution ):
	"""! @brief generate coverage histogram """
	
	values = []
	blocks = [ cov_values[ i : i + resolution ] for i in xrange( 0, len( cov_values ), resolution ) ]
	for block in blocks:
		values.append( min( [ np.mean( block ), saturation ] ) )
	
	fig, ax = plt.subplots()
	
	ax.hist( values, bins=300, color="lime" )
	ax.set_xlim( 0, 300 )
	
	ax.set_xlabel( "sequencing coverage depth" )
	ax.set_ylabel( "number of positions" )
	
	fig.savefig( outputfile, dpi=300 )
	plt.close( "all" )


def main( arguments ):
	"""! @brief runs everything """
	
	cov_file = arguments[ arguments.index( '--in' ) + 1 ]
	out_file = arguments[ arguments.index( '--out' ) + 1 ]
	
	if '--res' in arguments:
		resolution = int( arguments[ arguments.index( '--res' ) + 1 ] )
	else:
		resolution = 10000
	
	if '--sat' in arguments:
		saturation = int( arguments[ arguments.index( '--sat' ) + 1 ] )
	else:
		saturation = 300

	cov = load_cov( cov_file )
	
	# --- generate coverage histograms per chromosome --- #
	for key in sorted( cov.keys() )[:12]:
		outputfile = out_file + key + ".png"
		generate_hist( cov[ key ], outputfile, saturation, resolution )
	
	# --- generate per chromosome position coveage plot --- #
	generate_plot( cov, out_file, resolution, saturation )


if '--in' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
