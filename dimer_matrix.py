from Bio import SeqIO
import primer3
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
import math

n_primers = 0
all_primers = {}

for record in SeqIO.parse("primers_iniciales.fasta", "fasta"):
	primer_id = str(record.id)
	seq = str(record.seq).upper()
	all_primers[primer_id] = seq
	n_primers += 1 


names = all_primers.keys()

matrix = [[] for i in range((len(names)))]

for i in range(len(names)):
	for j in range(len(names)):
		matrix[i].append('hola')


df = pd.DataFrame(columns=names, index=names)
df_structures = pd.DataFrame(columns=names, index=names)

delta_g_total = 0

for i in names:
    for j in names:
        dimer = primer3.calcHeterodimer(all_primers[i], all_primers[j], output_structure=True)
        delta_g = dimer.dg/1000
        delta_g_total += delta_g
        if delta_g >= -5:
        	df.loc[i, j] = 0
        else:
        	df.loc[i, j] = delta_g

        hover = (	f'Primer 1: {i}' +
        			f'<br />Primer 2:{j}' +
        			f'<br />ΔG = {round(dimer.dg/1000, 2)}')
        try:
        	hover += (
        		f'<br />Estructura:'
        		f'<br />{dimer.ascii_structure_lines[0][4:]}'+
        		f'<br />{dimer.ascii_structure_lines[1][4:]}'+
        		f'<br />{dimer.ascii_structure_lines[2][4:]}'+
        		f'<br />{dimer.ascii_structure_lines[3][4:]}')
        except TypeError:
        	hover += f'<br />Sin estrctura'

        df_structures.loc[i, j] = hover


fig = go.Figure(
	data=go.Heatmap(
	z = df,
	zmax = 0,
	zmin = -20,
	text=df_structures,
	hoverinfo="text",
	hovertext=df_structures)
	)

fig.update_traces(hovertemplate=('<span style="font-size: 14px;'+
								'font-family: Courier">%{hovertext}</span>'))

	
fig.show()
print(delta_g_total)


#### PRIMERS FINALES

n_primers = 0
all_primers = {}

for record in SeqIO.parse("primers_finales.fasta", "fasta"):
	primer_id = str(record.id)
	seq = str(record.seq).upper()
	all_primers[primer_id] = seq
	n_primers += 1 


names = all_primers.keys()

matrix = [[] for i in range((len(names)))]

for i in range(len(names)):
	for j in range(len(names)):
		matrix[i].append('hola')


df = pd.DataFrame(columns=names, index=names)
df_structures = pd.DataFrame(columns=names, index=names)


delta_g_total = 0
for i in names:
    for j in names:
        dimer = primer3.calcHeterodimer(all_primers[i], all_primers[j], output_structure=True)
        delta_g = dimer.dg/1000
        delta_g_total += delta_g

        if delta_g >= -5:
        	df.loc[i, j] = 0
        else:
        	df.loc[i, j] = delta_g

        hover = (	f'Primer 1: {i}' +
        			f'<br />Primer 2:{j}' +
        			f'<br />ΔG = {round(dimer.dg/1000, 2)}')
        try:
        	hover += (
        		f'<br />Estructura:'
        		f'<br />{dimer.ascii_structure_lines[0][4:]}'+
        		f'<br />{dimer.ascii_structure_lines[1][4:]}'+
        		f'<br />{dimer.ascii_structure_lines[2][4:]}'+
        		f'<br />{dimer.ascii_structure_lines[3][4:]}')

        except TypeError:
        	hover += f'<br />Sin estrctura'

        df_structures.loc[i, j] = hover


fig = go.Figure(
	data=go.Heatmap(
	z = df,
	zmax = 0,
	zmin = -20,
	text=df_structures,
	hoverinfo="text",
	hovertext=df_structures)
	)

fig.update_traces(hovertemplate=('<span style="font-size: 14px;'+
								'font-family: Courier">%{hovertext}</span>'))

	
fig.show()
print(delta_g_total)