#!/usr/bin/env python3
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import sys
import os
import glob

class MetricMuncher:

	def __init__(self, paired_end=False):
		self.paired_end = paired_end
		self.busco = {}
		self.rnaquast = {}
		self.detonate = {}
		self.mapping = {}

	def read_transrate(self):
		"""
		TODO
		"""
		pass

	def read_busco(self, busco_file):
		assembly_name = busco_file.replace('.txt', '').replace('short_summary_', '')
		busco_metrics = {}
		with open(busco_file, 'r') as busco_reader:
			lines = busco_reader.readlines()
			busco_metrics['Complete BUSCOs'] = int(lines[9].split('\t')[1])
			busco_metrics['Complete and single-copy BUSCOs'] = int(lines[10].split('\t')[1])
			busco_metrics['Complete and duplicated BUSCOs'] = int(lines[11].split('\t')[1])
			busco_metrics['Fragmented BUSCOs'] = int(lines[12].split('\t')[1])
			busco_metrics['Missing BUSCOs'] = int(lines[13].split('\t')[1])
		self.busco[assembly_name] = busco_metrics

	def read_rnaquast(self, rnaquast_file):
		"""Reads the rnaquast short report. """
		rnaquast_metrics = pd.read_csv(rnaquast_file, sep='\t', index_col=0)
		print(rnaquast_metrics)
		# rnaquast_metrics = pd.to_numeric(rnaquast_file)

		# with open(rnaquast_file, 'r') as reader:
		# 	reader.readline()
		# 	self.rnaquast = self.extract_metrics(reader)
		# 	sensitivityFile = str(snakemake.input['rnaquast_sensitivity'])
		self.rnaquast = rnaquast_metrics

	@staticmethod
	def extract_metrics(metric_file_object, sep='\t'):
		"""Reads metrics from files such as rnaquast or detonate. Usually tab separated files."""
		metrics = {}
		for line in metric_file_object:
			metric, value = line.split(sep)
			metrics[metric] = float(value)
		return metrics
		
	def read_detonate(self, kc_file, contig_file, rsem_file):
		assembly_name = kc_file.replace('.txt', '').replace('kc_', '')
		detonate_metrics = {}
		with open(kc_file, 'r') as reader:
			detonate_metrics.update( self.extract_metrics(reader) )
		with open(rsem_file, 'r') as reader:
			detonate_metrics.update( self.extract_metrics(reader) )
		with open(contig_file, 'r') as reader:
			detonate_metrics.update( self.extract_metrics(reader) )
		self.detonate[assembly_name] = detonate_metrics
	
	def read_mapping(self, mapping_file):
		assembly_name = mapping_file.replace('_mapping_stats.txt', '')
		with open(mapping_file, 'r') as reader:
			percentage = float(reader.readline().rstrip('\n'))
			self.mapping[assembly_name] = {'Remapping rate': percentage}

	def create_all_metrics(self):
		"""
		Creates the final dataframe from all other metrics.
		"""

		detonate_df = pd.DataFrame.from_dict(self.detonate)
		mapping_df = pd.DataFrame.from_dict(self.mapping)
		busco_df = pd.DataFrame.from_dict(self.busco)
		frames = [detonate_df, mapping_df, self.rnaquast, busco_df]
		self.metrics = pd.concat(frames)

		# tools = ['DETONATE'] * len(detonate_df) + ['HISAT2'] * len(mapping_df) + ['rnaQUAST'] * len(self.rnaquast) + ['BUSCO'] * len(busco_df)
		# self.metrics.insert(0, 'Tool', tools)

	def save_all_metrics(self, file_name, sep='\t'):
		self.metrics.to_csv(file_name, sep=sep)

	def select_metrics(self):
		if self.paired_end:
			# rnaQUAST # TransRate # DETONATE # BUSCO
			metrics_to_extract = ['Database coverage', 'Duplication ratio', '95%-assembled isoforms', 'Misassemblies', 
								'mean_orf_percent', 'score', 'p_bases_uncovered',
								'kmer_compression_score', 'Score', 'unweighted_nucl_F1', 'unweighted_contig_F1',
								'Complete and single-copy BUSCOs', 'Complete and duplicated BUSCOs', 'Fragmented BUSCOs', 'Missing BUSCOs', 

								'Remapping rate']
								#								'Ex90N50', 
			self.high_values = ['Database coverage', '95%-assembled isoforms', 
							'mean_orf_percent', 'score',
							'Score', 'kmer_compression_score', 'unweighted_nucl_F1', 'unweighted_contig_F1', 
							'Complete and single-copy BUSCOs',

							'Remapping rate']
							#							'Ex90N50', 
			self.low_values = ['Duplication ratio', 'Misassemblies',
							'p_bases_uncovered',
							'Missing BUSCOs', 'Fragmented BUSCOs', 'Complete and duplicated BUSCOs']
		else:
			metrics_to_extract = ['Database coverage', 'Duplication ratio', '95%-assembled isoforms', 'Misassemblies', 
								'kmer_compression_score', 'Score', 'unweighted_nucl_F1', 'unweighted_contig_F1',
								'Complete and single-copy BUSCOs', 'Complete and duplicated BUSCOs', 'Fragmented BUSCOs', 'Missing BUSCOs', 

								'Remapping rate']
								#								'Ex90N50', 
			self.high_values = ['Database coverage', '95%-assembled isoforms', 
							'Score', 'kmer_compression_score', 'unweighted_nucl_F1', 'unweighted_contig_F1', 
							'Complete and single-copy BUSCOs',

							'Remapping rate']
							#							'Ex90N50', 
			self.low_values = ['Duplication ratio', 'Misassemblies',
							'Missing BUSCOs', 'Fragmented BUSCOs', 'Complete and duplicated BUSCOs']
		assert len(self.high_values) + len(self.low_values) == len(metrics_to_extract), 'Missing metrics in high/low arrays.'
		self.selected_metrics = self.metrics.loc[metrics_to_extract]

	def save_selected_metrics(self, file_name, sep='\t'):
		self.selected_metrics.to_csv(file_name, sep=sep)

	def normalize_selected_metrics(self):
		normalized_high_values = self.selected_metrics.loc[self.high_values].apply(normalize, axis='columns', raw=True, result_type='expand')
		normalized_low_values = self.selected_metrics.loc[self.low_values].apply(normalize, args=(True,), axis='columns', raw=True, result_type='expand')
		self.normalized_selected_metrics = normalized_high_values.append(normalized_low_values)

	def rename_selected_metrics(self):
		if self.paired_end:
			self.normalized_selected_metrics.rename(index={
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
			self.normalized_selected_metrics = self.normalized_selected_metrics.reindex(['Database\ncoverage', 'Duplication ratio', '95 %-assembled\nisoforms', 'Misassemblies', 
								'Mean ORF (%)', 'Assembly score', 'Uncovered bases (%)',
								'KC score', 'RSEM EVAL', 'Nucleotide F1', 'Contig F1',
								'BUSCOs (CS)', 'BUSCOs (CD)', 'BUSCOs (F)', 'BUSCOs (M)', 

								'Remapping rate'])
								#								'Ex90N50', 
		else:
			self.normalized_selected_metrics.rename(index={
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
			self.normalized_selected_metrics = self.normalized_selected_metrics.reindex(['Database\ncoverage', 'Duplication ratio', '95 %-assembled\nisoforms', 'Misassemblies',
							'KC score', 'RSEM EVAL', 'Nucleotide F1', 'Contig F1',
							'BUSCOs (CS)', 'BUSCOs (CD)', 'BUSCOs (F)', 'BUSCOs (M)', 

							'Remapping rate'])
							#							'Ex90N50', 

	def save_selected_normalized_metrics(self, file_name, sep='\t'):
		self.normalized_selected_metrics.to_csv(file_name, sep=sep)

	def calculate_column_sum(self):
		colSum = pd.DataFrame(self.normalized_selected_metrics.apply(sum, axis='rows')).transpose()
		colSum.index = ['SUM SCORE']
		self.final_metrics = pd.concat([self.normalized_selected_metrics, colSum])

	def save_heatmap(self, filename):
		'''Generates a heatmap for a pandas dataframe, but exclues the last row. For this script the last row is the sum score of all other rows and has another colour than the heatmap.
		
		Arguments:
			metrics_dataframe {DataFrame} -- dataframe to plot
			filename {str} -- file name of output file.
		'''
		sns.set()
		# Set the size of the plot to DIN A4 so values and axis are visible
		fig, ax = plt.subplots()
		# the size of A4 paper
		fig.set_size_inches(11.7, 8.27)

		# Set the mask to exclude the last row, which includes the sum score
		mask = np.zeros(self.final_metrics.shape)
		mask[-1,:] = True
		# Heatmap without the last row which contains the sum score
		hm = sns.heatmap(self.final_metrics,
						mask=mask, 
						cbar=True,
						cbar_kws={'pad': 0.01},
						annot=True,
						vmin=self.final_metrics.values[:-1,:].ravel().min(),
						vmax=self.final_metrics.values[:-1,:].ravel().max(),
						linewidth=0.75,
						cmap='coolwarm',
						fmt='.2f')

		# Set the mask to exlucde all data except of the sum score
		mask = np.ones(self.final_metrics.shape)
		mask[-1,:] = False
		# Heatmap for the last row.
		sns.heatmap(self.final_metrics,
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
		'''
		Normalizes any float/int array to the range of 0-1.
		:param array: float/int-type array
		:return: float-type array
		'''
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
	
	


mapping_stats = sorted(glob.glob('*_mapping_stats.txt'))
busco_stats = sorted(glob.glob('short_summary_*'))
rnaquast_stats = glob.glob('short_report.tsv')[0]
contig_stats = sorted(glob.glob('contig_nucl*'))
kc_stats = sorted(glob.glob('kc_*'))
rsem_stats = sorted(glob.glob('*.score'))


mm = MetricMuncher()
mm.read_rnaquast(rnaquast_stats)
for stats in mapping_stats:
	mm.read_mapping(stats)
for stats in busco_stats:
	mm.read_busco(stats)
for contig, kc, rsem in zip(contig_stats, kc_stats, rsem_stats):
	mm.read_detonate(kc, contig, rsem)


mm.create_all_metrics()
print(mm.metrics)
mm.save_all_metrics('all_metrics.csv')


mm.select_metrics()
mm.save_selected_metrics('selected_metrics.csv')
print(mm.selected_metrics)


mm.normalize_selected_metrics()
mm.rename_selected_metrics()
mm.save_selected_normalized_metrics('selected_normalized_metrics.csv')

mm.calculate_column_sum()
mm.save_heatmap('heatmap.svg')

