import argparse
import math
import random
import pandas as pd
import itertools
import re
import copy
import numpy as np
import plotly.express as px
import time
import Bio
from Bio import SeqIO
import subprocess
import incompatible_primers
import json
import pickle
import sys

testing_skip = False

pivot_flanking_size = 150

parser=argparse.ArgumentParser()


parser.add_argument("-in",
                    help="Input fasta file with the targets sequences.",
                    default = None,
                    type=str
                    )

parser.add_argument("-s0",
                    help=("Number of s0 sets (seeds) to generate." +
                        " Optimization will proceed from the best of these"),
                    default = 5,
                    type=int
                    )

parser.add_argument("-epochs",
                    help="Number of optimization cycles",
                    default = 10000,
                     type=int
                     )

parser.add_argument("-pools",
                    help="Number of pools",
                    default = 1,
                    type=int
                    )

parser.add_argument("-s",
                    help="max size of output amplicons",
                    default = 250,
                    type=int
                    )

parser.add_argument("-p",
                    help=("size of pivots regions to sequence."+
                        " Primers will be designed flanking this region"),
                    default = 100,
                    type=int
                    )

parser.add_argument("-pmin_size",
                    help="primers min size",
                    default = 18,
                    type=str
                    )

parser.add_argument("-pmax_size",
                    help="primers max size",
                    default = 25,
                    type=int
                    )

parser.add_argument("-pmin_gc",
                    help="primers min size",
                    default = .35,
                    type=float
                    )

parser.add_argument("-pmax_gc",
                    help="primers max size",
                    default = .7,
                    type=float
                    )

parser.add_argument("-verbose",
                    help="If verbose=1, the program will print the progress of the optimization",
                    default = 1,
                    type=int
                    )

parser.add_argument("-subseq",
                    help="Max size of the subsequences to be analyzed within the primers",
                    default = 8,
                    type=int
                    )

parser.add_argument("-salt",
                    help=(  "Salinity of the reaction medium. This parameter "+ 
                            "affects the hybridization ΔG of oligos."),
                    default = 0.18,
                    type=float
                    )

parser.add_argument("-temp",
                    help=(  "Temperature of hybridization in Celcius. This parameter "+ 
                            "affects the hybridization ΔG of oligos."),
                    default = 60,
                    type=float
                    )

parser.add_argument("-dg_ideal",
                    help="Optimal hybridization ΔG for oligos.",
                    default = -12.5,
                    type=float
                    )

parser.add_argument("-fchange",
                    help=(  "Initial fraction of primers pairs to be changed" +
                            " in the first optimization loop."),
                    default = .10,
                    type=float
                    )

parser.add_argument("-fdecay",
                    help=(  "Disminution per epoch of the fraction of primers"+ 
                            " pairs to be changed in each epoch"),
                    default = 0.0001,
                    type=float
                    )

parser.add_argument("-nmin",
                    help="Minimal amount of primers pairs to be changed in each epoch.",
                    default = 1,
                    type=int
                    )


args=vars(parser.parse_args())

# Archivos
targets_file = args['in']


# Parametros generales
num_s0 = args['s0']
epochs = args['epochs']
verbose = args['verbose']
num_pools = args['pools']

amplicon_max_size = args['s']
max_pivot_size = args['p']
num_pools = args['pools']

subseq_max_len = args['subseq']

# Parametros de los primers
primer_min_size = args['pmin_size']
primer_max_size = args['pmax_size']
primer_min_gc = args['pmin_gc']
primer_max_gc = args['pmax_gc']


# Parametros que afectan al ΔG de hibridacion
S = args['salt']
T = args['temp'] + 273.15
dg_ideal = args['dg_ideal']
primer_max_gc = args['pmax_gc']


# Parametros para un momentum gradient descent like optimization
frac_n_i = args['fchange']
frac_dis = args['fdecay']
n_min = args['nmin']


# ΔG para el inicio de la hibridizacion del oligo
dg_init = 2.09


# Carga y procesamiento de los datos termodinamicos utilizados

with open('data/thermo_parameters.json') as f:
    thermo_parameters = json.loads(f.read())

dg_dict = {}

for prop_seq in thermo_parameters.keys():
    dg = thermo_parameters[prop_seq]['dh'] 
    dg -= T * ((thermo_parameters[prop_seq]['ds'] + 0.368 * math.log(S))/1000)    
    dg_dict[prop_seq] = dg


def calc_primer_dg(primer_seq):

    dg = dg_init
    for i in range(len(primer_seq) - 1):
        dg += dg_dict[primer_seq[i:i+2]]

    return dg


def select_primers(seq_info_dict):

    first = random.choice(['fw','rv'])

    temp_forward_primers = copy.deepcopy(seq_info_dict['forward_primers'])
    temp_reverse_primers = copy.deepcopy(seq_info_dict['reverse_primers'])

    if first == 'fw':
        forward = random.choice(temp_forward_primers)

        while forward['primer_id'] in incompatible_set:
            # Si el forward seleccionado esta dentro del set de incompat.,
            # se remueve de la lista temporal de primers a elegir.
            
            temp_forward_primers.remove(forward)

            if len(temp_forward_primers) > 0:
                forward = random.choice(temp_forward_primers)
                # Si la lista se queda sin primers
                # osea, (len(temp_forward_primers) == 0)
                # se dejara el primer que ya estaba
            else:
                forward = seq_info_dict['forward']
                break


        reverse = random.choice(temp_reverse_primers)

        while (
                reverse['primer_id'] in incompatible_set or
                reverse['primer_5'] - forward['primer_5'] > amplicon_max_size
        ):

            temp_reverse_primers.remove(reverse)


            if len(temp_reverse_primers) > 0:
                reverse = random.choice(temp_reverse_primers)

            else:
                reverse = seq_info_dict['reverse']
                break

    elif first == 'rv':

        reverse = random.choice(temp_reverse_primers)

        while reverse['primer_id'] in incompatible_set :
            # Si el forward seleccionado esta dentro del set de incompat.,
            # se remueve de la lista temporal de primers a elegir.

            temp_reverse_primers.remove(reverse)

            if len(temp_reverse_primers) > 0:
                reverse = random.choice(temp_reverse_primers)
                # Si la lista se queda sin primers
                # osea, (len(temp_forward_primers) == 0)
                # se dejara el primer que ya estaba
            else:
                reverse = seq_info_dict['reverse']
                break



        forward = random.choice(temp_forward_primers)

        while (
                forward['primer_id'] in incompatible_set or
                reverse['primer_5'] - forward['primer_5'] > amplicon_max_size
        ):

            temp_forward_primers.remove(forward)

            if len(temp_forward_primers) > 0:
                forward = random.choice(temp_forward_primers)

            else:
                forward = seq_info_dict['forward']
                break

    seq_info_dict['forward'] = forward
    seq_info_dict['reverse'] = reverse


def gen_incompatible_set(pools_list):

    incompatible_set = set()

    for pool in pools_list:

        for seq_info in pool:
            try:
                forward_id = seq_info['forward']['primer_id']
                for element in incompatiblity_dict[forward_id]:
                    incompatible_set.add(element)
            except KeyError:
                pass

            try:
                reverse_id = seq_info['forward']['primer_id']
                for element in incompatiblity_dict[reverse_id]:
                    incompatible_set.add(element)
            except KeyError:
                pass

    return incompatible_set


def reverse_complement(seq):

    rev_comp_seq = []

    for base in seq[::-1]:
        if base == 'A':
            rev_comp_seq.append('T')
        elif base == 'T':
            rev_comp_seq.append('A')
        elif base == 'C':
            rev_comp_seq.append('G')
        elif base == 'G':
            rev_comp_seq.append('C')

    rev_comp_seq = ''.join(rev_comp_seq)

    return rev_comp_seq

def gen_subseq(largo):
    bases = ['A', 'C', 'T', 'G']
    yield from itertools.product(*([bases] * (largo+1)))


def distribute_list(seqs_info_list, num_pools):
    random.shuffle(seqs_info_list)  # shuffle the list randomly
    result = [[] for _ in range(num_pools)]  # create n empty sublists
    for i, item in enumerate(seqs_info_list):
        result[i % num_pools].append(item)  # distribute the items into the sublists in a round-robin fashion
    return result


"""## Creacion de diccionario de lista de diccionarios de secuencias seqs_info"""

seqs = []

for record in SeqIO.parse(targets_file, "fasta"):
  seq = str(record.seq).upper()
  id = str(record.id)

  # Extrayendo coordenadas de las secuencias
  coord_start, coord_end = id.split("-")
  coord_start = int(coord_start.split(":")[1])
  coord_end = int(coord_end)


  if len(seq) - pivot_flanking_size*2 > max_pivot_size :
    chunks = math.ceil((len(seq) - pivot_flanking_size*2)/max_pivot_size)
    for j in range(chunks):
      chunk_5 = max_pivot_size * j
      chunk_3 = (j+1)*max_pivot_size + pivot_flanking_size *2
      chunk_seq = seq[chunk_5:chunk_3]

      coord_start_chunk = coord_start + chunk_5
      coord_end_chunk = coord_start + chunk_3
      chunk_id = id + '_' +str(j+1)
      seqs.append((chunk_id, chunk_seq, coord_start_chunk, coord_end_chunk))

  else:
    seqs.append((id ,seq, coord_start, coord_end))

print(str(len(seqs)) + ' target subsequences generated')


"""Se genera una lista (seqs_info) donde cada elemento es un diccionario para almacenar la infomacion de cada secuencia (longitud, sitio de inicio y final del pivot).

Ademas, se le agrega a cada secuencia los sitos 5' maximos para los primers forward y reverse que se generaran.
Por ejemplo, para el sitio máximo 5' del primer forward (max_fw_5), se tiene en cuenta que cuando tenemos un reverse de tamaño minimo (primer_min_size) hibridando justo sobre el 3' del pivot, el primer forward solo podrá llegar a un sitio de manera que al calcular la distancia entre ambos primers sea el tamaño maximo del amplicon (amplicon_max_size). De esta manera evitamos que el programa genere primers que nunca se podrían utilizar por generar siempre amplicones de tamaños mayores al limite.
"""

seqs_info = []
for seq in seqs:
    seq_len = len(seq[1])
    pivot_start = pivot_flanking_size
    pivot_end = seq_len-pivot_flanking_size
    max_fw_5 =  pivot_end + primer_min_size - amplicon_max_size

    if max_fw_5 < 0:
      max_fw_5 = 0

    # Limite 5' para el primer fw para que usando un reverso de tamano minimo
    # no supere el maximo tamano de amplicon

    max_rv_5 = pivot_start - primer_min_size + amplicon_max_size
    if max_rv_5 > seq_len:
      max_rv_5 = seq_len

    # Limite 5' para el primer rv para que usando un forward de tamano minimo
    # no supere el maximo tamano de amplicon

    seqs_info.append({
    'seq' : seq[1], 'pivot_start' : pivot_start,
    'pivot_end': pivot_end, 'seq_len' : seq_len,
    'max_fw_5': max_fw_5, 'max_rv_5': max_rv_5,
    'seq_id': seq[0], 'coord_start' : seq[2],
    'coord_end' : seq[3]})


### Generacion de primers forward


for seq_info in seqs_info:

    '''
    generamos proto-primers. Estos seria una libreria de todas las secuencias
    desde el 3' del pivot hasta el 5' maximo, cada uno terminando en un 5'
    distinto.
    '''

    n = 1

    proto_primers = []
    seq = seq_info['seq']
    proto_primers_3 = seq_info['pivot_start']
    primer_5 = seq_info ['max_fw_5']
    primer_size = proto_primers_3 - primer_5

    while primer_size >= primer_min_size:
        primer = seq[primer_5:proto_primers_3]
        proto_primers.append({  'proto_primer_seq': primer,
                                'primer_5': primer_5})
        primer_5 += 1
        primer_size = proto_primers_3 - primer_5

    '''
    Ahora se editaran los proto-primers para que cumplan con los requisitos
    termodinamicos. Los optimo es un ΔG de hibridacion de aprox -11,5 kcal/mol.
    (aceptable de -10.5 a -12.5 kcal/mol).
    Para lograr que los primers tengan este ΔG se los va recortando a los
    proto-primers de a un nucleotido desde el 3' hasta que logran un ΔG < -12.5.
    Luego se sigue eliminando de a un nucleotido. Si la eliminacion hace que
    el ΔG se acerque aun mas al ΔG ideal (-11.5 kcal/mol), la eliminacion se
    acepta. Si la eliminacion aleja al ΔG del primer del ΔG ideal con respecto
    a antes de la eliminacion, se descarta la eliminacion.
    Una vez aceptado el primer, se guarda la secuencia del mismo y la posicion
    de su 5' en el diccionario de cada secuencia, dentro de la lista de
    secuencias seqs_info
    '''

    seq_info['forward_primers'] = []

    for primer in proto_primers:

        bases_recortadas = 0

        for i in range(len(primer['proto_primer_seq'])):

            dg_i = calc_primer_dg(primer['proto_primer_seq'])

            if dg_i < -12.5:
                primer['proto_primer_seq'] = primer['proto_primer_seq'][:-1]
                bases_recortadas += 1
                continue

            else:
                dist_i = abs(dg_i - dg_ideal)
                dg_f = calc_primer_dg(primer['proto_primer_seq'][:-1])
                dist_f = abs(dg_f - dg_ideal)

                if dist_f < dist_i:
                    primer['proto_primer_seq'] = primer['proto_primer_seq'][:-1]
                    bases_recortadas += 1

                else:
                    primer_seq = primer['proto_primer_seq']

                primer_gc = primer_seq.count('C') + primer_seq.count('G')
                primer_gc /= len(primer_seq)

                if primer_min_gc <= primer_gc <= primer_max_gc:

                  primer_id = seq_info['seq_id'] + '_primer_f_' + str(n)
                  n += 1
                  primer_3 = proto_primers_3 - bases_recortadas

                  primer_5_coord = seq_info['coord_start'] + primer['primer_5']
                  primer_3_coord = seq_info['coord_start'] + primer_3

                  seq_info['forward_primers'].append({'primer_seq': primer_seq,
                                                      'primer_5': primer['primer_5'],
                                                      'primer_3': primer_3,
                                                      'primer_id': primer_id,
                                                      'primer_5_coord' : primer_5_coord,
                                                      'primer_3_coord' : primer_3_coord})

                break

"""### Generacion de primers reverse

Se utiliza un sistema identico al de la generacion de primers forward (buscando en la otra direccion), con la diferencia que una vez generados los primers, estos se almacenan usando la funcion reverse_complement, que devuelve la secuencia reversa complementaria.

"""

### print('Generando primers reverse...')


for seq_info in seqs_info:

    n = 1

    proto_primers = []
    # set initial list for appending the proto-primers

    seq = seq_info['seq']
    proto_primers_3 = seq_info['pivot_end']
    primer_5 = seq_info['max_rv_5']
    primer_size = primer_5 - proto_primers_3


    while primer_size >= primer_min_size:
        primer = seq[proto_primers_3:primer_5]
        proto_primers.append({  'proto_primer_seq': primer,
                                'primer_5': primer_5})
        primer_5 -= 1
        primer_size = primer_5 - proto_primers_3

    seq_info['reverse_primers'] = []


    for primer in proto_primers:

        bases_recortadas = 0

        for i in range(len(primer['proto_primer_seq'])):

            dg_i = calc_primer_dg(primer['proto_primer_seq'])

            if dg_i < -12.5:
                primer['proto_primer_seq'] = primer['proto_primer_seq'][1:]
                bases_recortadas += 1
                continue

            else:
                dist_i = abs(dg_i - dg_ideal)
                dg_f = calc_primer_dg(primer['proto_primer_seq'][1:])
                dist_f = abs(dg_f - dg_ideal)

                if dist_f < dist_i:
                    primer_seq = reverse_complement(primer['proto_primer_seq'][1:])
                    bases_recortadas += 1

                else:
                    primer_seq = reverse_complement(primer['proto_primer_seq'])

                primer_gc = primer_seq.count('C') + primer_seq.count('G')
                primer_gc /= len(primer_seq)

                if primer_min_gc <= primer_gc <= primer_max_gc:
                  primer_id = seq_info['seq_id'] + '_primer_r_' + str(n)
                  n += 1
                  primer_3 = proto_primers_3 + bases_recortadas

                  primer_5_coord = seq_info['coord_start'] + primer['primer_5']
                  primer_3_coord = seq_info['coord_start'] + primer_3

                  seq_info['reverse_primers'].append({'primer_seq': primer_seq,
                                                      'primer_5': primer['primer_5'],
                                                      'primer_3': primer_3,
                                                      'primer_id': primer_id,
                                                      'primer_5_coord' : primer_5_coord,
                                                      'primer_3_coord' : primer_3_coord})


                break


"""### Filtro de primers sin pareja posible

Para evitar que al momento de seleccionar los primers se seleccione un primer que no tenga un par correspondiente (fw/rv) cuya combinacion genere un amplicon de un tamano permitido (< max_amplicon_size), buscamos los extremos 5' de los fw y rv mas cercanos a los pivots. Todos aquellos primers cuyos 5' tengan una distancia mayor al max_amplicon_size con respecto a estos extremos, seran filtrados.

Si esto no se incluye, al momento de seleccionar los primers, si se toma un primer sin pareja posible, la funcion queda atrapada en el while loop que busca a un primer compatible.
"""

# print('Filtrando primers sin pareja posible')


for seq_info in seqs_info:

  '''Buscar el 5' del primer rv mas cercano al pivot
  Seteamos el 5' mas cercano, al extremo (max_rv_5) para q vaya bajando
  '''
  closest_rv_5 = seq_info['max_rv_5']

  for primer in seq_info['reverse_primers']:
    if primer['primer_5'] < closest_rv_5:
      closest_rv_5 = primer['primer_5']

  '''Buscar el 5' del primer fw mas cercano al pivot

  Seteamos el 5' mas cercano, al extremo (max_rv_5) para q vaya bajando
  '''

  closest_fw_5 = seq_info['max_fw_5']

  for primer in seq_info['forward_primers']:
    if primer['primer_5'] > closest_fw_5:
      closest_fw_5 = primer['primer_5']

  '''Eliminamos de la lista de fw aquellos cuya diferencia con el
  5' del rv mas ceracano al pivot sea mayor al amplicon_max_size.
  Iteramos la lista en sentido inverso (reversed), para no generar
  conflictos en la iteracion al eliminar elementos
  '''
  for primer in reversed(seq_info['forward_primers']):
    if closest_rv_5 - primer['primer_5'] > amplicon_max_size:
      seq_info['forward_primers'].remove(primer)

  for primer in reversed(seq_info['reverse_primers']):
    if primer['primer_5'] - closest_fw_5 > amplicon_max_size:
      seq_info['reverse_primers'].remove(primer)


# print('Chequeando la existencia de targets sin primers')

"""Chequeo que no exista ningun target del cual no se hayan obtenido primers"""

no_forwards = []
to_pop = []

for i, seq_info in enumerate(seqs_info):
  if seq_info['forward_primers'] == []:
    # print(seq_info['seq_id'], 'sin primers!')
    no_forwards.append(seq_info)
    to_pop.append(i)

no_reverse = []

for i, seq_info in enumerate(seqs_info):
  if seq_info['reverse_primers'] == []:
    # print(seq_info['seq_id'], 'sin primers!')
    no_reverse.append(seq_info)
    to_pop.append(i)


for index in sorted(to_pop, reverse=True):
    del seqs_info[index]

if len(to_pop) > 0:
  print(len(to_pop), 'target subsequences without possible primers!')
else:
  print('All target subsequences have possible primers!')

######################## GENERACION DE TABLAS DE INCOMPATIBILDIAD

with open('all_primers.fasta', 'w') as file:
  for seq_info in seqs_info:
    for rev_primer in seq_info['reverse_primers']:
        file.write('>' + rev_primer['primer_id'] + '\n')
        file.write(rev_primer['primer_seq'] + '\n')

    for fw_primer in seq_info['forward_primers']:
        file.write('>' + fw_primer['primer_id'] + '\n')
        file.write(fw_primer['primer_seq'] + '\n')


####### Mapeo de primers al genoma de referencia

    mapping = False

    if mapping:

        print('Mapeando primers al genoma de referencia...')

        bwa_aln_args = ['bwa', 'aln',
                        '-n', '1',
                        '-N',
                        '-t', '2',
                        '-M', '1',
                        'example/hg19.fasta',
                        'all_primers.fasta']


        with open('example/results_aln.sai', 'w') as fp:
            bwa_aln_cmd = subprocess.run(bwa_aln_args, stdout=fp)

        bwa_samse_args = ['bwa','samse',
                        '-n', '10000000',
                        '-f', 'out.sam',
                        'example/hg19.fasta',
                        'example/results_aln.sai',
                        'all_primers.fasta']


        bwa_samse_cmd = subprocess.run(bwa_samse_args)



    ### se remueve el header del sam file (lineas iniciadas con @)
        grep_args = ['grep', '-v', '@', 'out.sam']

        with open('out-mod.sam', 'w') as fp:
            grep_cmd = subprocess.run(grep_args, stdout=fp)

        print('Alineamiento terminado')

        incompatible_primers.create_df('out-mod.sam')
        incompatiblity_dict = incompatible_primers.find_incompatibilities()
        with open('testing/incompatiblity_dict.txt', 'w') as convert_file:
             convert_file.write(json.dumps(incompatiblity_dict))
    else:
        with open('testing/incompatiblity_dict.txt') as json_file:
            incompatiblity_dict = json.load(json_file)


"""### Seleccion aleatoria de primers.

Se utiliza la funcion select_primers, la cual toma un par de primers aleatorio para cada target y los almacena el el diccionario de cada secuencia target con las keys 'forward' y 'reverse'
"""

'''
Calculo del Loss del set inicial L0

Creamos una lista con todas las posibles subsecuencias (hashs) que se puedan
formar con el numero de bases máximo setado anteriormente (subseq_max_len).
Por ejemplo, si subseq_max_len = 8, la lista estara compuesta por:

['A', 'T', 'G', 'C', 'AA', 'AT', 'AC', ... , 'GGGGGGGC', 'GGGGGGGG']
'''

sub_seqs = []

for i in range(subseq_max_len):
  for x in gen_subseq(i):
      sub_seqs.append(''.join(x))

## distribuimos las seq_info en pooles aleatoriamente



pools_list = distribute_list(seqs_info, num_pools)

S0_list = []

pools_hash = [0 for _ in range(num_pools)]

for s in range(num_s0):

    print('Seed S0 #' + str(s+1))

    total_loss = 0

    for p_index, pool in enumerate(pools_list):
        all_primers = []
        for seq_info in pool:
            # Se crea el incompatible_set vacio para evtiar errores al selecccionar
            # los primers
            incompatible_set = set()

            select_primers(seq_info)

            '''Generamos una lista (all_primers) del set de primers inicial seleccionado
            aleatoriamente S0'''

            # print(seq_info)
            all_primers.append(seq_info['forward']['primer_seq'])
            all_primers.append(seq_info['reverse']['primer_seq'])


            '''Creamos la hashtable utilizando las subsecuencias generadas en el paso
            anterior como hashs. Inicializamos todos los valores de cada hash con 0.

            Luego, se hace buscan estas subsecuencias en todos los primers. En cada match,
            se suma al valor del hash de cada subsecuencia H(subseq) el valor 1/(1+d), donde
            d es la distancia de la subsecuencia la 3' del primer.

            Se implementa un loop while de manera de encontrar todas las veces que un hash
            esta en un primer. De manera contraria solo se computaria el primer match y se
            ignorarian los siguientes. '''



        seqs_hash = {x: 0 for x in sub_seqs}

        for key in seqs_hash.keys():
            for primer in all_primers:
                res = primer.find(key)
                while res > -1:
                    dist_3 = len(primer) - res - len(key)
                    seqs_hash[key] += 1/(1+dist_3)

                    # Searches for other appearing of the hash in the same primer

                    primer = primer[res+1:]
                    res = primer.find(key)

        '''En el siguiente paso se debe buscar la presencia del reverso complementario
        de cada subsecuencia en los primers. Para esto generamos una lista con el
        reverso complementario de todos los primers de nuestro set S0.
        '''

        rev_comp_primers = []
        for primer in all_primers:
            rev_comp_primers.append(reverse_complement(primer))

        '''Cuando se encuentra un match, se debe calcula H(subsecuencia)/(1/d), donde d
        es la distancia al 3' del reverso complementario.

        Para esto comenzamos primero generamos un diccionario vacio que tiene como keys
        todo los k_mers (subsecuencias de largo 1 a subseq_max_len) y a cada una de esta
        s se les asigna como valor una lista vacia.
        '''

        k_mers_hash = {}

        for i in range(subseq_max_len):
            for primer in rev_comp_primers:
                for n in range(len(primer)-i):
                    k_mer = primer[n:i+n+1]
                    k_mers_hash[k_mer] = []

        '''Luego barremos todas las subsecuencias de los primers nuevamente, agregando a
        la lista de cada k_mer (encontrada k_mers_hash[k_mer], las distancias al 3' del
        primer de cada subsecuencia.'''

        for i in range(subseq_max_len):
            for primer in rev_comp_primers:
                for n in range(len(primer)-i):
                    k_mer = primer[n:i+n+1]
                    dist_3 = len(primer)-(n+i+1)
                    k_mers_hash[k_mer].append(dist_3)

        '''
        Luego, para cada subsecuencia calculamos el badness con la siguiente ecuacion:


        $H(s) = sum_{}^{}\frac{H(rev.comp)}{d+1} \times 2^{len}\times 2^{numGC}$

        Cada H(s) es adicionado al Loss, obteniendo el resultado final al final del
        bucle.'''

        pool_loss = 0

        for k_mer, distances in k_mers_hash.items():
            k_mer_len = len(k_mer)
            gc = k_mer.count('C') + k_mer.count('G')
            H = seqs_hash[k_mer]
            pool_loss += sum([H/(dist_3+1) for dist_3 in distances]) * 2**k_mer_len * 2**gc

            pools_hash[p_index] = seqs_hash

        # print('Pool', p_index+1, '- L(S0)=',  round(pool_loss))

        total_loss += pool_loss

    print('L(S0) = '+ str(round(total_loss)))


    S0_list.append((total_loss, copy.deepcopy(pools_list), copy.deepcopy(pools_hash)))


"""Se selecciona como set S0 a usarse como seed, a aquel que presente menor valor de L(S0)"""

S0_losses = []

for total_loss, _, _ in S0_list:
  S0_losses.append(total_loss)

index_min_loss = np.argmin(S0_losses)

total_loss = S0_list[index_min_loss][0]

print('Lowest L(S0) = ' + str(round(total_loss)))

pools_list = S0_list[index_min_loss][1]
pools_hash = S0_list[index_min_loss][2]

pools_list_S0 = S0_list[index_min_loss][1]

# Generamos un record de los total_loss para luego graficarlos.

lossS0 = total_loss
loss_record = []
loss_record.append(total_loss)


"""## Bucle de optimizacion

Para optimizar el set de primers, en cada generacion **g** se genera temporalmente un set S<sub>g</sub> cambiando alazar un par de primers de alguno de los targets.
Luego, se debe recalcular el hashtable restando el aporte de los primers de la generacion anterior g-1 y agregando los nuevos primers seleccionados.
Finalmente se utiliza esta nueva hash table recaluclada para hacer el calculo del loss L(S<sub>g</sub>), de la misma manera que se calculo L(S<sub>0</sub>).

Si L(S<sub>g</sub>) > L(S<sub>g-1</sub>), se descarta el set de primers S<sub>g</sub> y se vuelve a establecer el set S<sub>g-1</sub>

Si L(S<sub>g</sub>) < L(S<sub>g-1</sub>), el set de primers S<sub>g</sub> pasa a ser el nuevo set.

"""

'''Indices de las secuencias para luego hacer el sampleo aleatorio para la
optimizacion
'''

lowest_loss = total_loss

try:

    for m in range(epochs):

        pools_sizes = [len(pool) for pool in pools_list]

        pools_indexes = []
        for pool in pools_list:
            indexes = [n for n in range(len(pool))]
            pools_indexes.append(indexes)


        total_loss_i = total_loss

        start = time.time()


        pools_list_temp = copy.deepcopy(pools_list)
        pools_hash_temp = copy.deepcopy(pools_hash)

        end = time.time()

        # print('tiempo deepcopy', end - start)


        ##### Seleccion random de un par de primers a cambiar

        seq_totales = sum([len(pool) for pool in pools_list_temp])
        n_change = (frac_n_i-m*frac_dis)*seq_totales
        n_change = round(n_change)


        if n_change < n_min:
          n_change = n_min

        if verbose >= 2:
            print('Secuencias a cambiar = ' + str(n_change) + '/' +  str(seq_totales))


         ##### Generando lista de incompatibilidades

        incompatible_set = gen_incompatible_set(pools_list_temp)

        ########### Sampleo aleatorio de indices de secuencias a modif #############
        
        # Se seleccionan los pooles a cambiar, dandole un peso segun el tamano de cada pool

        pool_to_change = random.choices(range(len(pools_list_temp)), weights=pools_sizes, k=n_change)

        # se crea la lista de random_indexes, donde se guardaran los indices de
        # las secuencias a cambiar. Por ejemplo, si hay 3 pooles, podemos tener
        # random_indexes = [[1, 4, 7, 21], [4, 10, 8], [22, 1]], lo que quiere
        # decir que se cambiaran las seuencias 1, 4, 7 y 21 del primer pool (pool 0),
        # las secuencias 4, 10 y 8 del pool 1, etc.

        random_indexes = [[] for pool in pools_list_temp]


        for pool_index in pool_to_change:
            try:
                seq_index = random.choice(pools_indexes[pool_index])
                random_indexes[pool_index].append(seq_index)
                pools_indexes[pool_index].remove(seq_index)
            except IndexError as e:
                print(e)
                print(pool_index)
                print(pools_indexes)

        for pool_random_indexes in random_indexes:
            pool_random_indexes.sort(reverse=True)


        ############# Modificando la hashtable #############

        #### Restar los H correspondientes a los k_mers de estos primers eliminados

        # lista para alamacenar los nuevos primers en cada pool,
        # para luego sumar su contribucion a sus respectivos
        # hashtables

        pools_new_primers = [[] for pool in pools_list_temp]

        for p_index, pool in enumerate(pools_list_temp):

            del_primers = []

            # Dentro del pool p_index, se seleccionan las secuencias listadas
            # en random_indexes[p_index] (seleccionadas para eliminarse de ese
            # pool) y se agregan los primers correspondinetes a estas secuencias
            # a una lista de del_primers a eliminarse.

            for random_index in random_indexes[p_index]:
              del_forward = pool[random_index]['forward']['primer_seq']
              del_primers.append(del_forward)
              del_reverse = pool[random_index]['reverse']['primer_seq']
              del_primers.append(del_reverse)


            for key in pools_hash_temp[p_index].keys():
                for primer in del_primers:
                    res = primer.find(key)
                    while res > -1:
                        dist_3 = len(primer) - res - len(key)
                        pools_hash_temp[p_index][key] -= 1/(1+dist_3)
                        if pools_hash_temp[p_index][key] < 1E-3:
                            pools_hash_temp[p_index][key] = 0
                        primer = primer[res+1:]
                        res = primer.find(key)

            

            for i, random_index in enumerate(random_indexes[p_index]):

                ## Se asignan nuevos primers a las secuencias en las cuales
                ## se eliminaron sus primers en el paso anterior.
                
                select_primers(pool[random_index])

                assigned_pool = random.randint(0, len(pools_list_temp)-1)

                if assigned_pool != p_index:
                    seq_to_move = pool.pop(random_index)
                    pools_list_temp[assigned_pool].append(seq_to_move)
                else:
                    seq_to_move = pool[random_index]

                new_forward = seq_to_move['forward']['primer_seq']
                pools_new_primers[assigned_pool].append(new_forward)
                new_reverse = seq_to_move['reverse']['primer_seq']
                pools_new_primers[assigned_pool].append(new_reverse)

                # asignar los pares de primers aleatoriamente a alguno de los pooles


                #### Sumar los H correspondientes a los k_mers de los nuevos primers

        total_loss = 0

        for p_index, pool in enumerate(pools_list_temp):

            new_primers = pools_new_primers[p_index]

            for key in pools_hash_temp[p_index].keys():
                for primer in new_primers:
                    res = primer.find(key)
                    while res > -1:
                        dist_3 = len(primer) - res - len(key)
                        pools_hash_temp[p_index][key] += 1/(1+dist_3)
                        primer = primer[res+1:]
                        res = primer.find(key)


            ##### Regenerando lista de primers

            all_primers_temp = []

            for seq_info in pool:
                all_primers_temp.append(seq_info['forward']['primer_seq'])
                all_primers_temp.append(seq_info['reverse']['primer_seq'])



            ######### Calculando el Loss de el nuevo set L(Sg) #########


            rev_comp_primers = []
            for primer in all_primers_temp:
                rev_comp_primers.append(reverse_complement(primer))


            k_mers_hash = {}

            for i in range(subseq_max_len):
                for primer in rev_comp_primers:
                    for n in range(len(primer)-i):
                        k_mer = primer[n:i+n+1]
                        k_mers_hash[k_mer] = []

            for i in range(subseq_max_len):
                for primer in rev_comp_primers:
                    for n in range(len(primer)-i):
                        k_mer = primer[n:i+n+1]
                        dist_3 = len(primer)-(n+i+1)
                        k_mers_hash[k_mer].append(dist_3)

            end = time.time()
            pool_loss = 0

            for k_mer, distances in k_mers_hash.items():
                k_mer_len = len(k_mer)
                gc = k_mer.count('C') + k_mer.count('G')
                H = pools_hash_temp[p_index][k_mer]
                pool_loss += sum([H/(dist_3+1) for dist_3 in distances]) * 2**k_mer_len * 2**gc

                if verbose >= 2:
                    print('Total loss - L(S)=', round(total_loss))
                    print('Lowest loss - L(S)=', round(lowest_loss))

            if verbose >= 2:
                print('Pool', p_index+1 ,'- L(S)=', round(pool_loss))

            total_loss += pool_loss

        if verbose >= 1:

            epoch_msg = ('Epoch = ' + str(m+1) + '/' + str(epochs) + '\n'
                        'Loss = ' + str(round(total_loss)) + '\n'
                        'Lowest Loss = ' + str(round(lowest_loss))+ '\n')

            print(epoch_msg, end='\033[3A')



        loss_record.append(total_loss)


        ''' Calcular p, la probabilidad de aceptar un set que empeora el loss.
        Esto nos permitiria evitar caer en minimos locales

        accept_S_temp tomara el valor True o False con una probabilidad dependiente
        de la p'''

        prob_accept = False

        if prob_accept == True:
          C0 = lossS0/500
          k = .1
          C = C0 * math.exp(-m*k)
          try:
            p = min(1, math.exp((total_loss_i-loss)/(C*num_seq)))
          except OverflowError:
            p = 1
          finally:
            # Si (Lg) < (Lg-1) => p = 1 y siempre se aceptara el cambio
            accept_S_temp = random.random() < p

        else:
          if total_loss_i > total_loss:
            accept_S_temp = True
          else:
            accept_S_temp = False


        if accept_S_temp:
          pools_hash = copy.deepcopy(pools_hash_temp)
          pools_list = copy.deepcopy(pools_list_temp)
          lowest_loss = total_loss
          lowest_loss_epoch = m+1

          count_incompat = True

          all_primers_id = []
          for pool in pools_list_temp:
            for seq_info in pool:
                try:
                    all_primers_id.append(seq_info['forward']['primer_id'])
                    all_primers_id.append(seq_info['reverse']['primer_id'])
                except:
                    print(seq_info)
                    raise TypeError

          if count_incompat:
              inc_c = 0
              for primer_id in all_primers_id:
                  try:
                      for incomp_primer in incompatiblity_dict[primer_id]:
                          if incomp_primer in incompatible_set:
                              inc_c += 1
                  except KeyError:
                      pass
          
          if verbose >= 2:
            print('Incompatibilidades =', inc_c)

        else:
          total_loss = total_loss_i

except KeyboardInterrupt:
    print("Optimization finalizada por el usuario")

finally:

    ##### Guardado de resultados ####

    with open('testing/seqs_info_inicial.pkl', 'wb') as file:
        pickle.dump(pools_list_S0, file)

    with open('testing/seqs_info_final.pkl', 'wb') as file:
        pickle.dump(pools_list, file)

    with open('testing/loss_record.pkl', 'wb') as file:
        pickle.dump(loss_record, file)


    """## Grafico del Loss vs Epochs"""

    """## Set de primers iniciales"""


    with open('primers_iniciales.fasta', 'w+') as file:

      for index, pool in enumerate(pools_list_S0):
        file.write(('------------'+' Primers Pool ' + str(index + 1)+ ' ------------'))
        file.write('\n')
        for seq_info in pool:
          file.write(('>' + seq_info['seq_id'] + 'fw' + '\n' + seq_info['forward']['primer_seq']))
          file.write('\n')
          file.write(('>' + seq_info['seq_id'] + 'rv' + '\n' + seq_info['reverse']['primer_seq']))
          file.write('\n')


    with open('primers_finales.fasta', 'w+') as file:

      for index, pool in enumerate(pools_list):
        file.write(('------------'+' Primers Pool ' + str(index + 1)+ ' ------------'))
        file.write('\n')
        for seq_info in pool:
            file.write(('>' + seq_info['seq_id'] + 'fw' + '\n' + seq_info['forward']['primer_seq']))
            file.write('\n')
            file.write(('>' + seq_info['seq_id'] + 'rv' + '\n' + seq_info['reverse']['primer_seq']))
            file.write('\n')


    with open('primers_finales.bed', 'w+') as file:
        
        header = 'track name=Primers description="Primers Designed by Moo" itemRgb="On"'
        file.write((header))
        file.write('\n')

        RGB_colors =     [  '255,105,97', 
                            '255,180,128',
                            '248,243,141',
                            '66,214,164',
                            '8,202,209',
                            '89,173,246',
                            '157,148,255',
                            '199,128,232']

        rgb_cicle = itertools.cycle(RGB_colors)


        for pool in pools_list:
            for seq_info, rgb_code in zip(pool, rgb_cicle):

                fw_5 = str(seq_info['forward']['primer_5_coord'])
                rv_5 = str(seq_info['reverse']['primer_5_coord'])
                fw_size = str(len(seq_info['forward']['primer_seq']))
                rv_size = str(len(seq_info['reverse']['primer_seq']))

                rv_3 = str(seq_info['reverse']['primer_3_coord'])
                relative_rv_3 = str(int(rv_3) - int(fw_5))


                seq_id = seq_info['seq_id']
                line = ('chr7' + ' ' 
                        + fw_5 + ' '
                        + rv_5 + ' '
                        + seq_id + ' '
                        + '0' + ' '
                        + '+' + ' '
                        + fw_5 + ' '
                        + rv_5 + ' '
                        + rgb_code + ' '
                        + '2' + ' '
                        + fw_size + ','
                        + rv_size + ' '
                        + '0,'
                        + relative_rv_3
                        )

                file.write(line)
                file.write('\n')
