## **Lipidomics Project**

### Getting Started
R 3.6.0 or higher version and RStudio needed.

### Goal
The project is designed as a referece for lipidomics experiments. All the data are collected via Lipid Search software. 
The pipeline will implement quality control, statistics analysis and visualization for the mass spectrometry based data.

### Introduction
<ul>
<li>FWL_lipidomics2.8.R</li>
<li>FWL_lipidomics.functions2.8.R</li>
<li>FWL_FattyAcids.functions.R</li>
<li>Archive (old versions of pipeline)
	<ul>...</ul>
</li>
<li>convert_data</li>
<li>fattyAcids_saturation_analysis
	<ul>
	<li>README.md</li>
	<li>fattyAcids_saturation_analysis2.1.r</li>
	<li>fattyAcids_saturation_analysis.FUNCTIONS2.1.r</li>
	</ul>
</li>
<li>test_data
	<ul>
		<li>sample1.csv</li>
		<li>sample1_raw.csv</li>
	</ul>
</li>
<li>pathway_visualization
	<ul>
		<li>README.md</li>
		<li>Def_Graphic_Map.txt</li>
		<li>LiPa_Shiny_V1.0_g.R</li>
		<li>LiPa_code_V0.8.R</li>
		<li>Standards_to_use_pmol.txt</li>
		<li>cleandata.csv</li>
		<li>lipid_map.png</li>
	</ul>
</li>
<li>quality_control
	<ul>
		<li>README.md</li>
		<li>script.R</li>
	</ul>
</li>
<li>statistics_quantification
	<ul>
		<li>Heatmap_dendrogram_plot.R</li>
		<li>LipidSearch_statistics_two_factors.R</li>
		<li>Volcano.R</li>
	</ul>
</li>

The main pipeline name starts with FWL and it will implement a consistent analysis for the search files, including fatty acids saturation analysis. And the saturation analysis could also be implemented independently via pipeline which has loose data filtering standard in fattyAcids_saturation_analysis folder. The main pipeline is built on the previous work of qualitiy_control and statistics_quantification files. 
