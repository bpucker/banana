### Boas Pucker ###
### bpucker@cebitec.uni-bielefeld.de ###
### v0.2 ###

### adjusted for banana based on https://doi.org/10.1371/journal.pone.0164321 ###

__usage__ = """
					python SnpEff_result_parser.py
					--in <FULL_PATH_TO_INPUT_VCF>
					--out <FULL_PATH_TO_OUTPUT_TEXT_FILE>
					--anno <ANNOTATION_FILE>
					
					bug reports and feature requests: bpucker@cebitec.uni-bielefeld.de
					"""

import re, sys

# ---- end of imports --- #


def find_high_impact_variants( input_file, output_file, annotation ):
	"""! @brief load only first high impact effect per gene """
	
	ended_genes = []
	
	effect_impact_positions = { }
	#collect effect positions to avoid two high impact effects at same SNP
	
	splice_variants_counter = 0
	premature_stop_counter = 0
	frame_shift_counter = 0
	lost_stop_counter = 0
	
	high_counter = 0
	moderate_counter = 0
	low_counter = 0
	modifier_counter = 0
	
	premature_stops = []
	
	data_for_extraction = []
	
	total_annotated_variant_counter = 0
	total_small_variant_counter = 0
	total_snp_counter = 0
	
	with open( input_file, "r" ) as f:
		line = f.readline()
		while line:
			if line[0] != '#':
				parts = line.strip().split('\t')
				total_annotated_variant_counter += 1
				try:
					x = abs( len( parts[3] ) - len( parts[4] ) )
				except IndexError:
					print line
				if abs( len( parts[3] ) - len( parts[4] ) ) < 100:	#filter out large InDels
					total_small_variant_counter += 1
					if len( parts[3] ) == len( parts[4] ):
						total_snp_counter += 1
					subparts = parts[7].split(',')
					for subpart in subparts:
						# --- check all variant annotations and select only one (highest impact and first) --- #
						try:
							ID = re.findall( "Ma[0-9Un_random]+_g\d+", subpart )[0]
							if ID not in ended_genes:
								# --- analyse high impact effects --- #
								if 'HIGH' in line:
									try:
										position = effect_impact_positions[ parts[0]+parts[1] ]									
									except:
										effect_annotation = "ERROR"
										if 'stop_gained' in subpart:
											premature_stop_counter += 1
											premature_stops.append( ID )
											ended_genes.append( ID )
											effect_annotation = 'stop_gained'
										elif 'stop_lost' in subpart:
											lost_stop_counter += 1
											effect_annotation = 'stop_lost'
										elif 'splice_region_variant' in subpart:
											splice_variants_counter += 1
											effect_annotation = 'splice_region_variant'
										elif 'frameshift' in subpart:
											frame_shift_counter += 1
											effect_annotation = 'frameshift'
										
										effect_impact_positions.update( { parts[0]+parts[1]: "" } )
										high_counter += 1
										if effect_annotation != "ERROR":
											data_for_extraction.append( [ parts[0], parts[1], parts[3], parts[4], re.findall( "Ma[0-9Un_random]+_g\d+", ID )[0], effect_annotation ] )
										
								elif 'MODERATE' in line:
									try:
										position = effect_impact_positions[ parts[0]+parts[1] ]											
									except:
										effect_impact_positions.update( { parts[0]+parts[1]: "" } )
										moderate_counter += 1
								elif 'LOW' in line:
									try:
										position = effect_impact_positions[ parts[0]+parts[1] ]											
									except:
										effect_impact_positions.update( { parts[0]+parts[1]: "" } )
										low_counter += 1
								elif 'MODIFIER' in line:
									try:
										position = effect_impact_positions[ parts[0]+parts[1] ]											
									except:
										effect_impact_positions.update( { parts[0]+parts[1]: "" })
										modifier_counter += 1
						except:
							pass	#print "ERROR: " + subpart
			line = f.readline()
	print "RESULTS:"
	
	
	print "total annotated variants: " + str( total_annotated_variant_counter )
	print "total small variants: " + str( total_small_variant_counter )
	print "total number of SNPs: " + str( total_snp_counter )
	
	
	print "number of splice variants: " + str( splice_variants_counter )
	print "number of premature stops: " + str( premature_stop_counter )
	print "number of frame shifts: " + str( frame_shift_counter )
	print "number of lost stops: " + str( lost_stop_counter )
	
	
	
	print "HIGH: " + str( high_counter )
	print "MODERATE: " + str( moderate_counter )
	print "LOW: " + str( low_counter )
	print "MODIFIER: " + str( modifier_counter )
	
	print "premature stop check: " + str( len( premature_stops ) )
	print "premature stop check (unique): " + str( len( list( set( premature_stops ) ) ) )
	
	
	lengths_of_ref_alleles = []
	lengths_of_alt_alleles = []
	# --- write collected data into output file --- #
	with open( output_file, "w" ) as out:
		header = [ "Chromosome", "Position", "ReferenceAllel", "AlternativeAllel", "GeneID", "EffectType", "Annotation" ]
		out.write( "\t".join( header ) + '\n' )
		for each in data_for_extraction:
			try:
				out.write( "\t".join( each + [ annotation[ each[4] ] ] )  + '\n' )
			except KeyError:
				out.write( "\t".join( each ) + '\tn/a\n' )
			lengths_of_ref_alleles.append( len( each[2] ) )
			lengths_of_alt_alleles.append( len( each[3] ) )
	
	#print sorted( lengths_of_ref_alleles )[-20:]
	#print max( lengths_of_alt_alleles )


def load_annotation( annotation_file ):
	"""! @brief load functional gene annotation from given file """
	
	annotation = {}
	with open( annotation_file, "r" ) as f:
		line = f.readline()
		while line:
			parts = line.strip().split('\t')
			annotation.update( { parts[0]: "_%_".join( parts[1:] ) } )
			line = f.readline()
	return annotation


def main( arguments ):
	"""! @brief run everything """
	
	input_file = arguments[ arguments.index('--in')+1 ]
	output_file = arguments[ arguments.index('--out')+1 ]
	if '--anno' in arguments:
		annotation_file = arguments[ arguments.index('--anno')+1 ]
		annotation = load_annotation( annotation_file )
	else:
		annotation = {}
	
	find_high_impact_variants( input_file, output_file, annotation )

if __name__ == '__main__':
	
	if '--in' in sys.argv and '--out' in sys.argv:
		main( sys.argv )
	else:
		sys.exit( __usage__ )
