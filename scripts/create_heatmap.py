import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def generate_heatmap(metrics_dataframe, filename):
	"""Generates a heatmap for a pandas dataframe, but exclues the last row. For this script the last row is the sum score of all other rows and has another colour than the heatmap.
	
	Arguments:
		metrics_dataframe {DataFrame} -- dataframe to plot
		filename {str} -- file name of output file.
	"""
	sns.set()
	# Set the size of the plot to DIN A4 so values and axis are visible
	fig, ax = plt.subplots()
	# the size of A4 paper
	fig.set_size_inches(11.7, 8.27)

	# Set the mask to exclude the last row, which includes the sum score
	mask = np.zeros(metrics_dataframe.shape)
	mask[-1,:] = True
	# Heatmap without the last row which contains the sum score
	hm = sns.heatmap(metrics_dataframe,
					mask=mask, 
					cbar=True,
					cbar_kws={'pad': 0.01},
					annot=True,
					vmin=metrics_dataframe.values[:-1,:].ravel().min(),
					vmax=metrics_dataframe.values[:-1,:].ravel().max(),
					linewidth=0.75,
					cmap='coolwarm',
					fmt='.2f')

	# Set the mask to exlucde all data except of the sum score
	mask = np.ones(metrics_dataframe.shape)
	mask[-1,:] = False
	# Heatmap for the last row.
	sns.heatmap(metrics_dataframe,
					mask=mask,
					alpha=0.5,
					cbar=False,
					linewidth=0.75,
					cmap = 'binary',
					vmax=0,
					annot=True,
					annot_kws={'color':'black', 'weight': 'bold'},
					fmt='.2f')

	plt.xticks(rotation=45, ha='right')		# Rotate the xaxis at 45 degree
	plt.tight_layout()                      # command that x and y labels are not cut off.
	plt.savefig(filename, format='svg')		# Save the figure as svg

def normalize(array, reverse=False):
	"""
	Normalizes any float/int array to the range of 0-1.
	:param array: float/int-type array
	:return: float-type array
	"""
	if reverse:
		minimum = max(array)
		maximum = min(array)
	else:
		maximum = max(array)
		minimum = min(array)
	returnArray = []
	for a in array:
		if (maximum - minimum) == 0:
			returnArray.append(0)
		else:
			returnArray.append(abs((a-minimum)/(maximum-minimum)))
	return(returnArray)

assembler = snakemake.params['assemblies'] # Used assemblers as given in the contig file

# Create empty dataframe with assemblers as columns
assemblerDict = {}
for a in assembler:
	assemblerDict[a] = []
infoDF = pd.DataFrame(data=assemblerDict)

##############################
########## rnaQUAST ##########
##############################
sensitivityFile = str(snakemake.input['rnaquast_sensitivity'])

# whitespace delimiter, but only if its 2 or more whitespaces after another
sensitivity = pd.read_csv(sensitivityFile, delimiter=r"\s{2,}", index_col=0, na_values = '*', engine='python')
sensitivity = sensitivity.dropna(axis=0)
infoDF = pd.concat([infoDF, sensitivity])

misassembliesFile = str(snakemake.input['rnaquast_misassemblies'])
misassemblies = pd.read_csv(misassembliesFile, delimiter=r"\s{2,}", index_col=0, na_values = '*', engine='python')
misassemblies = misassemblies.dropna(axis=0)
infoDF = pd.concat([infoDF, misassemblies])

##############################
########## Transrate #########
##############################
paired_end = True
try:
	assemblies = str(snakemake.input['transrate'])
except:
	paired_end = False
if paired_end:
	transrateInfos = pd.read_csv(assemblies, index_col=0).transpose()
	transrateInfos = transrateInfos.dropna(axis=0)
	transrateInfos.columns = assembler
	infoDF = pd.concat([infoDF, transrateInfos])

##############################
########## DETONATE ##########
##############################
detonateDir = str(snakemake.input['detonateDir'])
# kc_paths = snakemake.input['kc']
# rsem_eval_paths = snakemake.input['rsem_eval']

detonateInfos = pd.DataFrame()

for count, a in enumerate(assembler):
	#kc_path = kc_paths[count]
	kc_path = detonateDir + '/kc_' + a + '.txt'
	kc = pd.read_csv(kc_path, delimiter='\t', index_col=0, names=[a])

	contig_nucl_path = detonateDir + '/contig_nucl_' + a + '.txt'
	contig = pd.read_csv(contig_nucl_path, delimiter='\t', index_col=0, names=[a])
	all = kc.append(contig)

	#rsem_eval_path = rsem_eval_paths[count]
	rsem_eval_path = detonateDir + '/rsem_eval_' + a + '.score'
	rsem = pd.read_csv(rsem_eval_path, delimiter='\t', index_col=0, names=[a])
	all = all.append(rsem)

	detonateInfos = pd.concat([detonateInfos, all], axis=1)

infoDF = pd.concat([infoDF, detonateInfos])

##############################
########### BUSCO ############
##############################
busco_dir = str(snakemake.input['busco'])
buscoInfos = pd.DataFrame()
for count, a in enumerate(assembler):
	busco_results = f'{busco_dir}/run_{a}/short_summary_{a}.txt'
	with open(busco_results, 'r') as busco_reader:
		lines = busco_reader.readlines()
		complete_buscos = int(lines[9].split('\t')[1])
		complete_single_buscos = int(lines[10].split('\t')[1])
		complete_duplicated_buscos = int(lines[11].split('\t')[1])
		fragemented_buscos = int(lines[12].split('\t')[1])
		missing_buscos = int(lines[13].split('\t')[1])
	metric_names = ['Complete BUSCOs', 'Complete and single-copy BUSCOs', 'Complete and duplicated BUSCOs', 'Fragmented BUSCOs', 'Missing BUSCOs']
	metrics = [complete_buscos, complete_single_buscos, complete_duplicated_buscos, fragemented_buscos, missing_buscos]
	col_name = [a]
	temp_busco_results = pd.DataFrame(metrics, columns=col_name, index=metric_names)
	buscoInfos = pd.concat([buscoInfos, temp_busco_results], axis=1)

infoDF = pd.concat([infoDF, buscoInfos])

##############################
########### Ex90N50 ##########
##############################
EX90N50_dir = str(snakemake.input['EX90N50'])
EX90N50_infos = pd.DataFrame()
for count, a in enumerate(assembler):
	stats_file = f'{EX90N50_dir}/{a}/ExN50.stats'
	with open(stats_file, 'r') as reader:
		for line in reader:
			if line.startswith('90'):
				EX90N50_value = int(line.split()[2])
	metric_names = ['Ex90N50']
	metrics = [EX90N50_value]
	col_name = [a]
	temp_results = pd.DataFrame(metrics, columns=col_name, index=metric_names)
	EX90N50_infos = pd.concat([EX90N50_infos, temp_results], axis=1)

infoDF = pd.concat([infoDF, EX90N50_infos])

##############################
####### Remapping rates ######
##############################

mapping_rates = str(snakemake.input['mapping_rates'])
mapping_infos = pd.read_csv(mapping_rates, sep = '\t', index_col = 0)

infoDF = pd.concat([infoDF, mapping_infos])

# Remove rows that contain NaN values
infoDF = infoDF.dropna(axis=0)

# Remove rnaspades
infoDF.drop(columns = ['rnaSPAdes_cl'], inplace=True)
# Replace names in headers
infoDF.rename(columns={'karma_p': 'karma*', 'Grouper_p': "Grouper*", "combined": "Combined"}, inplace=True)
##############################
######### SAVE FILES #########
##############################

# Save unnormalized
rawFile = str(snakemake.output.raw)

infoDF.to_csv(rawFile, sep='\t')


# Save chosen values in etra file.
reducedRawFile = str(snakemake.output.minRaw)
if paired_end:
	# rnaQUAST # TransRate # DETONATE # BUSCO
	metricsToExtract = ['Database coverage', 'Duplication ratio', '95%-assembled isoforms', 'Misassemblies', 
						'mean_orf_percent', 'score', 'p_bases_uncovered',
						'kmer_compression_score', 'Score', 'unweighted_nucl_F1', 'unweighted_contig_F1',
						'Complete and single-copy BUSCOs', 'Complete and duplicated BUSCOs', 'Fragmented BUSCOs', 'Missing BUSCOs', 
						'Ex90N50', 
						'mapping_rate_per_base']
	highValues = ['Database coverage', '95%-assembled isoforms', 
					'mean_orf_percent', 'score',
					'Score', 'kmer_compression_score', 'unweighted_nucl_F1', 'unweighted_contig_F1', 
					'Complete and single-copy BUSCOs',
					'Ex90N50', 
					'mapping_rate_per_base']
	lowValues = ['Duplication ratio', 'Misassemblies',
					'p_bases_uncovered',
					'Missing BUSCOs', 'Fragmented BUSCOs', 'Complete and duplicated BUSCOs']
else:
	metricsToExtract = ['Database coverage', 'Duplication ratio', '95%-assembled isoforms', 'Misassemblies', 
						'kmer_compression_score', 'Score', 'unweighted_nucl_F1', 'unweighted_contig_F1',
						'Complete and single-copy BUSCOs', 'Complete and duplicated BUSCOs', 'Fragmented BUSCOs', 'Missing BUSCOs', 
						'Ex90N50', 
						'mapping_rate_per_base']
	highValues = ['Database coverage', '95%-assembled isoforms', 
					'Score', 'kmer_compression_score', 'unweighted_nucl_F1', 'unweighted_contig_F1', 
					'Complete and single-copy BUSCOs',
					'Ex90N50', 
					'mapping_rate_per_base']
	lowValues = ['Duplication ratio', 'Misassemblies',
					'Missing BUSCOs', 'Fragmented BUSCOs', 'Complete and duplicated BUSCOs']
assert len(highValues) + len(lowValues) == len(metricsToExtract), 'Missing metrics in high/low arrays.'

reducedRawScores = infoDF.loc[metricsToExtract]
reducedRawScores.to_csv(reducedRawFile, '\t')


normDF_highValues = reducedRawScores.loc[highValues].apply(normalize, axis='columns', raw=True, result_type='expand')
normDF_lowValues = reducedRawScores.loc[lowValues].apply(normalize, args=(True,), axis='columns', raw=True, result_type='expand')
normDF = normDF_highValues.append(normDF_lowValues)

# Rename some indices
if paired_end:
	normDF.rename(index={
				'Database coverage': 'Database\ncoverage',
				'Score': "RSEM EVAL",
                'kmer_compression_score': 'KC score',
				'mean_orf_percent': 'Mean ORF (%)',
				'score': 'Assembly score',
				'p_bases_uncovered': 'Uncovered bases (%)',
				'unweighted_nucl_F1': 'Nucleotide F1',
				'unweighted_contig_F1': 'Contig F1',
				'mapping_rate_per_base': 'Remapping rate',
				'95%-assembled isoforms': '95 %-assembled\nisoforms',
				'Complete and single-copy BUSCOs': 'BUSCOs (CS)',
				'Complete and duplicated BUSCOs': 'BUSCOs (CD)',
				'Fragmented BUSCOs': 'BUSCOs (F)',
				'Missing BUSCOs': 'BUSCOs (M)'
				}, inplace=True)
	normDF = normDF.reindex(['Database\ncoverage', 'Duplication ratio', '95 %-assembled\nisoforms', 'Misassemblies', 
						'Mean ORF (%)', 'Assembly score', 'Uncovered bases (%)',
						'KC score', 'RSEM EVAL', 'Nucleotide F1', 'Contig F1',
						'BUSCOs (CS)', 'BUSCOs (CD)', 'BUSCOs (F)', 'BUSCOs (M)', 
						'Ex90N50', 
						'Remapping rate'])
else:
	normDF.rename(index={
			'Database coverage': 'Database\ncoverage',
			'Score': "RSEM EVAL",
			'kmer_compression_score': 'KC score',
			'p_bases_uncovered': 'Uncovered bases (%)',
			'unweighted_nucl_F1': 'Nucleotide F1',
			'unweighted_contig_F1': 'Contig F1',
			'mapping_rate_per_base': 'Remapping rate',
			'95%-assembled isoforms': '95 %-assembled\nisoforms',
			'Complete and single-copy BUSCOs': 'BUSCOs (CS)',
			'Complete and duplicated BUSCOs': 'BUSCOs (CD)',
			'Fragmented BUSCOs': 'BUSCOs (F)',
			'Missing BUSCOs': 'BUSCOs (M)'
			}, inplace=True)
	normDF = normDF.reindex(['Database\ncoverage', 'Duplication ratio', '95 %-assembled\nisoforms', 'Misassemblies',
					'KC score', 'RSEM EVAL', 'Nucleotide F1', 'Contig F1',
					'BUSCOs (CS)', 'BUSCOs (CD)', 'BUSCOs (F)', 'BUSCOs (M)', 
					'Ex90N50', 
					'Remapping rate'])



colSum = pd.DataFrame(normDF.apply(sum, axis='rows')).transpose()
colSum.index = ['SUM SCORE']
normDF = pd.concat([normDF, colSum])

# Save reduced normalized
normalizedFile = str(snakemake.output.normalized)
normDF.columns = infoDF.columns
normDF.to_csv(normalizedFile, sep='\t')
heatmapSVG = str(snakemake.output.heatmap)
generate_heatmap(normDF, filename=heatmapSVG)



# Save all normalized
# normalizedFile = str(snakemake.output.normalized)
# normDF = infoDF.apply(normalize, axis='columns', raw=True, result_type='expand')
# normDF.columns = assembler
# normDF.to_csv(normalizedFile, sep='\t')
