import math
import random
import pandas as pd
import itertools
import re
import copy
import numpy as np
import plotly.express as px
import time
import subprocess
import argparse
import Bio
from Bio import SeqIO

argParser = argparse.ArgumentParser()
argParser.add_argument("-t", "--targets", help="fasta file with targets sequences")

args = argParser.parse_args()
print("args=%s" % args)

print("args.name=%s" % args.targets)

targets_fasta = args.targets

intron_extra = 150

"""## Parametros

### Parametros de optimizacion

*  num_s0 = numero de sets S0 que se generan. La optimizacion comenzara del mejor de estos.
*  epochs = numero de ciclos de optimizacion
*  verbose = si se hace un print del avance (Si es True, se muestran los prints)
"""

num_s0 = 5
epochs = 5000
verbose = True

"""Determinacion de parmetros del amplicon y el pivot (region a secuenciar)"""

amplicon_max_size = 300
max_pivot_size = 100
primer_min_size = 18
primer_max_size = 30

"""Parametros para la generacion de secuencias random para el testing"""

num_seq = 48
min_gc = 35
max_gc = 70
min_len = 200
max_len = 400

"""Parametros para el calculo del Badness.
El parámetro seq_max_len es el tamaño máximo de subsecuencia que se analizará
dentro de los primers, buscando complementaridad (dimerizacion) con los
otros primers

"""

seq_max_len = 8

"""Parametros termodinamicos para caluclar ΔG hibridacion primers. Se obtuvieron
 los valors de ΔH y ΔS de hibridacion para cada dinucleotido (Santalucia J Jr
 et al., 2004)"""

thermo_parameters = {'AA' : {'dh' : -7.6, 'ds' : -21.3},
                     'AT' : {'dh' : -7.2, 'ds' : -20.4},
                     'AC' : {'dh' : -8.4, 'ds' : -22.4},
                     'AG' : {'dh' : -7.8, 'ds' : -21},
                     'TA' : {'dh' : -7.2, 'ds' : -21.3},
                     'TT' : {'dh' : -7.6, 'ds' : -21.3},
                     'TC' : {'dh' : -8.2, 'ds' : -22.2},
                     'TG' : {'dh' : -8.5, 'ds' : -22.7},
                     'CA' : {'dh' : -8.5, 'ds' : -22.7},
                     'CT' : {'dh' : -7.8, 'ds' : -21},
                     'CC' : {'dh' : -8.0, 'ds' : -19.9},
                     'CG' : {'dh' : -10.6, 'ds' : -27.2},
                     'GA' : {'dh' : -8.2, 'ds' : -22.2},
                     'GT' : {'dh' : -8.4, 'ds' : -22.4},
                     'GC' : {'dh' : -9.8, 'ds' : -24.4},
                     'GG' : {'dh' : -8.0, 'ds' : -19.9}
                     }

"""Algunos parametros que afectan al ΔG"""

S = 0.18 #salinity
T = 60 + 273.15 # temp in K
dg_init = 2.09
dg_ideal = -12.5

primer_min_gc = .35
primer_max_gc = .70

"""Calculo de ΔG de hibridacion afectados por salindad (S) y temperatura (T)"""

dg_dict = {}

for prop_seq in thermo_parameters.keys():
    dg = thermo_parameters[prop_seq]['dh'] - T * ((thermo_parameters[prop_seq]['ds'] + 0.368 * math.log(S))/1000)
    dg_dict[prop_seq] = dg

def random_seq(length, gc):

    bases = ['A', 'T', 'C', 'G']
    weights = [(1-gc)/2, (1-gc)/2, gc/2, gc/2]

    random_bases = random.choices(
                        population = bases,
                        weights=weights,
                        k=length)

    return ''.join(random_bases)


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


def calc_primer_dg(primer_seq):

    dg = dg_init
    for i in range(len(primer_seq) - 1):
        dg += dg_dict[primer_seq[i:i+2]]

    return dg


def select_primers(seq_info_dict):

    first = random.choice(['fw','rv'])
    c = 0

    if first == 'fw':
        forward = random.choice(seq_info_dict['forward_primers'])
        reverse = random.choice(seq_info_dict['reverse_primers'])

        while reverse['primer_5'] - forward['primer_5'] > amplicon_max_size:
            reverse = random.choice(seq_info_dict['reverse_primers'])

    elif first == 'rv':
        reverse = random.choice(seq_info_dict['reverse_primers'])
        forward = random.choice(seq_info_dict['forward_primers'])

        while reverse['primer_5'] - forward['primer_5'] > amplicon_max_size:
            forward = random.choice(seq_info_dict['forward_primers'])

    seq_info_dict['forward'] = forward['primer_seq']
    seq_info_dict['reverse'] = reverse['primer_seq']


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

"""## Creacion de diccionario de lista de diccionarios de secuencias seqs_info"""

seqs = []

for record in SeqIO.parse(targets_fasta, "fasta"):
  seq = str(record.seq).upper()
  id = str(record.id)
  print
  if len(seq) - intron_extra*2 > max_pivot_size :
    chunks = math.ceil(len(seq)/max_pivot_size)
    for j in range(chunks):
      chunk_5 =j*max_pivot_size - intron_extra
      if chunk_5 < 0:
        chunk_5 = 0
      else:
        pass
      chunk_3 = (j+1)*max_pivot_size + intron_extra
      chunk_seq = seq[chunk_5:chunk_3]
      chunk_id = id + '_' +str(j+1)
      seqs.append((chunk_id, chunk_seq))

  else:
    seqs.append((id ,seq))

"""Se genera una lista (seqs_info) donde cada elemento es un diccionario para almacenar la infomacion de cada secuencia (longitud, sitio de inicio y final del pivot).

Ademas, se le agrega a cada secuencia los sitos 5' maximos para los primers forward y reverse que se generaran.
Por ejemplo, para el sitio máximo 5' del primer forward (max_fw_5), se tiene en cuenta que cuando tenemos un reverse de tamaño minimo (primer_min_size) hibridando justo sobre el 3' del pivot, el primer forward solo podrá llegar a un sitio de manera que al calcular la distancia entre ambos primers sea el tamaño maximo del amplicon (amplicon_max_size). De esta manera evitamos que el programa genere primers que nunca se podrían utilizar por generar siempre amplicones de tamaños mayores al limite.
"""

seqs_info = []
for seq in seqs:
    seq_len = len(seq[1])
    pivot_start = intron_extra
    pivot_end = seq_len-intron_extra

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
    'seq_id': seq[0]
    })

"""## Generacion de primers

### Generacion de primers forward
"""

for seq_info in seqs_info:

    '''
    generamos proto-primers. Estos seria una libreria de todas las secuencias
    desde el 3' del pivot hasta el 5' maximo, cada uno terminando en un 5'
    distinto.
    '''

    proto_primers = []
    seq = seq_info['seq']
    primers_3 = seq_info['pivot_start']
    primer_5 = seq_info ['max_fw_5']
    primer_size = primers_3 - primer_5

    while primer_size >= primer_min_size:
        primer = seq[primer_5:primers_3]
        proto_primers.append({'proto_primer_seq': primer, 'primer_5': primer_5})
        primer_5 += 1
        primer_size = primers_3 - primer_5

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

        for i in range(len(primer['proto_primer_seq'])):

            dg_i = calc_primer_dg(primer['proto_primer_seq'])

            if dg_i < -12.5:
                primer['proto_primer_seq'] = primer['proto_primer_seq'][:-1]
                continue

            else:
                dist_i = abs(dg_i - dg_ideal)
                dg_f = calc_primer_dg(primer['proto_primer_seq'][:-1])
                dist_f = abs(dg_f - dg_ideal)

                if dist_f < dist_i:
                    primer['proto_primer_seq'] = primer['proto_primer_seq'][:-1]

                else:
                    primer_seq = primer['proto_primer_seq']

                primer_gc = primer_seq.count('C') + primer_seq.count('G')
                primer_gc /= len(primer_seq)
                if primer_min_gc <= primer_gc <= primer_max_gc:

                  seq_info['forward_primers'].append({'primer_seq': primer_seq,
                                                      'primer_5': primer['primer_5']})
                break

"""### Generacion de primers reverse

Se utiliza un sistema identico al de la generacion de primers forward (buscando en la otra direccion), con la diferencia que una vez generados los primers, estos se almacenan usando la funcion reverse_complement, que devuelve la secuencia reversa complementaria.
"""

for seq_info in seqs_info:

    proto_primers = []
    # set initial list for appending the proto-primers

    seq = seq_info['seq']
    primers_3 = seq_info['pivot_end']
    primer_5 = seq_info['max_rv_5']
    primer_size = primer_5 - primers_3


    while primer_size >= primer_min_size:
        primer = seq[primers_3:primer_5]
        proto_primers.append({'proto_primer_seq': primer, 'primer_5': primer_5})
        primer_5 -= 1
        primer_size = primer_5 - primers_3

    seq_info['reverse_primers'] = []


    for primer in proto_primers:

        for i in range(len(primer['proto_primer_seq'])):

            dg_i = calc_primer_dg(primer['proto_primer_seq'])

            if dg_i < -12.5:
                primer['proto_primer_seq'] = primer['proto_primer_seq'][1:]
                continue

            else:
                dist_i = abs(dg_i - dg_ideal)
                dg_f = calc_primer_dg(primer['proto_primer_seq'][1:])
                dist_f = abs(dg_f - dg_ideal)

                if dist_f < dist_i:
                    primer['proto_primer_seq'] = primer['proto_primer_seq'][1:]

                else:
                    primer_seq = reverse_complement(primer['proto_primer_seq'])

                primer_gc = primer_seq.count('C') + primer_seq.count('G')
                primer_gc /= len(primer_seq)
                if primer_min_gc <= primer_gc <= primer_max_gc:

                   seq_info['reverse_primers'].append({'primer_seq': primer_seq,
                                                      'primer_5': primer['primer_5']})
                break

"""### Filtro de primers sin pareja posible

Para evitar que al momento de seleccionar los primers se seleccione un primer que no tenga un par correspondiente (fw/rv) cuya combinacion genere un amplicon de un tamano permitido (< max_amplicon_size), buscamos los extremos 5' de los fw y rv mas cercanos a los pivots. Todos aquellos primers cuyos 5' tengan una distancia mayor al max_amplicon_size con respecto a estos extremos, seran filtrados.

Si esto no se incluye, al momento de seleccionar los primers, si se toma un primer sin pareja posible, la funcion queda atrapada en el while loop que busca a un primer compatible.
"""

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

"""Chequeo que no exista ningun target del cual no se hayan obtenido primers"""

no_forwards = []
for i, seq_info in enumerate(seqs_info):
  if seq_info['forward_primers'] == []:
    print(i)
    print(seq_info['seq_id'], 'sin primers!')
    no_forwards.append(seq_info)
    seqs_info.pop(i)

no_reverse = []
for i, seq_info in enumerate(seqs_info):
  if seq_info['reverse_primers'] == []:
    print(i)
    print(seq_info['seq_id'], 'sin primers!')
    no_reverse.append(seq_info)
    seqs_info.pop(i)

len(seqs_info)

"""### Seleccion aleatoria de primers.

Se utiliza la funcion select_primers, la cual toma un par de primers aleatorio para cada target y los almacena el el diccionario de cada secuencia target con las keys 'forward' y 'reverse'
"""

S0_list = []

for s in range(num_s0):
  for seq_info in seqs_info:
      select_primers(seq_info)
  '''Generamos una lista (all_primers) del set de primers inicial  seleccionado
  aleatoriamente S0'''

  all_primers = []
  for seq_info in seqs_info:
      all_primers.append(seq_info['forward'])
      all_primers.append(seq_info['reverse'])
  '''
  Calculo del Loss del set inicial L0

  Creamos una lista con todas las posibles subsecuencias (hashs) que se puedan
  formar con el numero de bases máximo setado anteriormente (seq_max_len).
  Por ejemplo, si seq_max_len = 8, la lista estara compuesta por:

  ['A', 'T', 'G', 'C', 'AA', 'AT', 'AC', ... , 'GGGGGGGC', 'GGGGGGGG']
  '''

  sub_seqs = []

  for i in range(seq_max_len):
      for x in gen_subseq(i):
          sub_seqs.append(''.join(x))
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
  todo los k_mers (subsecuencias de largo 1 a seq_max_len) y a cada una de esta
  s se les asigna como valor una lista vacia.
  '''

  k_mers_hash = {}

  for i in range(seq_max_len):
      for primer in rev_comp_primers:
          for n in range(len(primer)-i):
              k_mer = primer[n:i+n+1]
              k_mers_hash[k_mer] = []
  '''Luego barremos todas las subsecuencias de los primers nuevamente, agregando a
  la lista de cada k_mer (encontrada k_mers_hash[k_mer], las distancias al 3' del
  primer de cada subsecuencia.'''

  for i in range(seq_max_len):
      for primer in rev_comp_primers:
          for n in range(len(primer)-i):
              k_mer = primer[n:i+n+1]
              dist_3 = len(primer)-(n+i+1)
              k_mers_hash[k_mer].append(dist_3)
  '''
  Luego, para cada subsecuencia calculamos el badness con la siguiente ecuacion:


  $H(s) = \sum_{}^{}\frac{H(rev.comp)}{d+1} \times 2^{len}\times 2^{numGC}$

  Cada H(s) es adicionado al Loss, obteniendo el resultado final al final del
  bucle.'''

  loss = 0

  for k_mer, distances in k_mers_hash.items():
      k_mer_len = len(k_mer)
      gc = k_mer.count('C') + k_mer.count('G')
      H = seqs_hash[k_mer]
      loss += sum([H/(dist_3+1) for dist_3 in distances]) * 2**k_mer_len * 2**gc

  print('L(S0)=',  round(loss))

  S0_list.append((loss, copy.deepcopy(seqs_info), copy.deepcopy(seqs_hash)))

"""Se selecciona como set S0 a usarse como seed, a aquel que presente menor valor de L(S0)"""

S0_losses = []

for loss, seq_info, seqs_hash in S0_list:
  S0_losses.append(loss)

index_min_loss = np.argmin(S0_losses)

loss = S0_list[index_min_loss][0]
seqs_info = S0_list[index_min_loss][1]
seqs_hash = S0_list[index_min_loss][2]

# Generamos un record de los loss para luego graficarlos.

lossS0 = loss
loss_record = []
loss_record.append(loss)

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
indexes = [n for n in range(len(seqs_info))]

# Fraccion de secuencias a cambiar durante el primer ciclo
frac_n_i = .05
# Numero de secuencias a cambiar durante el primer ciclo

# Disminucion de la fraccion de secuencias en cada ciclo
frac_dis = 0.0001

# Minimo de secuencias a cambiar en cada ciclo
n_min = 1

for m in range(epochs):

    if verbose == True:
      print('Epoch = ' + str(m+1) + '/' + str(epochs) )

    loss_i = loss

    start = time.time()


    seqs_info_temp = copy.deepcopy(seqs_info)
    seqs_hash_temp = copy.deepcopy(seqs_hash)

    end = time.time()
    # print('tiempo deepcopy', end - start)


    ##### Seleccion random de un par de primers a cambiar

    n_change = (frac_n_i-m*frac_dis)*len(seqs_info_temp)
    n_change = round(n_change)

    if n_change < n_min:
      n_change = n_min

    # print('Cambiando', n_change, 'primers')

    random_indexes = random.sample(indexes, n_change)

    ############# Modificando la hashtable #############

    #### Restar los H correspondientes a los k_mers de estos primers eliminados
    start = time.time()

    del_primers = []

    for random_index in random_indexes:
      del_forward = seqs_info_temp[random_index]['forward']
      del_primers.append(del_forward)
      del_reverse = seqs_info_temp[random_index]['reverse']
      del_primers.append(del_reverse)


    for key in seqs_hash_temp.keys():
        for primer in del_primers:
            res = primer.find(key)
            while res > -1:
                dist_3 = len(primer) - res - len(key)
                seqs_hash_temp[key] -= 1/(1+dist_3)
                if seqs_hash_temp[key] < 1E-3:
                    seqs_hash_temp[key] = 0
                primer = primer[res+1:]
                res = primer.find(key)


    for random_index in random_indexes:
      select_primers(seqs_info_temp[random_index])

    end = time.time()
    # print('tiempo resta H primers eliminados', end - start)



    #### Sumar los H correspondientes a los k_mers de los nuevos primers
    start = time.time()

    new_primers = []

    for random_index in random_indexes:
      new_forward = seqs_info_temp[random_index]['forward']
      new_primers.append(new_forward)
      new_reverse = seqs_info_temp[random_index]['reverse']
      new_primers.append(new_reverse)


    for key in seqs_hash_temp.keys():
        for primer in new_primers:
            res = primer.find(key)
            while res > -1:
                dist_3 = len(primer) - res - len(key)
                seqs_hash_temp[key] += 1/(1+dist_3)
                primer = primer[res+1:]
                res = primer.find(key)

    end = time.time()
    # print('tiempo suma H primers nuevos', end - start)


    ##### Regenerando lista de primers
    all_primers_temp = []
    for seq_info in seqs_info_temp:
        all_primers_temp.append(seq_info['forward'])
        all_primers_temp.append(seq_info['reverse'])


    ######### Calculando el Loss de el nuevo set L(Sg) #########

    start = time.time()

    rev_comp_primers = []
    for primer in all_primers_temp:
        rev_comp_primers.append(reverse_complement(primer))


    k_mers_hash = {}

    for i in range(seq_max_len):
        for primer in rev_comp_primers:
            for n in range(len(primer)-i):
                k_mer = primer[n:i+n+1]
                k_mers_hash[k_mer] = []

    for i in range(seq_max_len):
        for primer in rev_comp_primers:
            for n in range(len(primer)-i):
                k_mer = primer[n:i+n+1]
                dist_3 = len(primer)-(n+i+1)
                k_mers_hash[k_mer].append(dist_3)

    end = time.time()
    # print('tiempo HASH rev comp', end - start)

    loss = 0

    for k_mer, distances in k_mers_hash.items():
        k_mer_len = len(k_mer)
        gc = k_mer.count('C') + k_mer.count('G')
        H = seqs_hash_temp[k_mer]
        loss += sum([H/(dist_3+1) for dist_3 in distances]) * 2**k_mer_len * 2**gc

    print('L(S)=', round(loss))
    loss_record.append(loss)

    end = time.time()
    # print('tiempo suma Loss', end - start)



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
        p = min(1, math.exp((loss_i-loss)/(C*num_seq)))
      except OverflowError:
        p = 1
      finally:
        # Si (Lg) < (Lg-1) => p = 1 y siempre se aceptara el cambio
        accept_S_temp = random.random() < p

    else:
      if loss_i > loss:
        accept_S_temp = True
      else:
        accept_S_temp = False


    if accept_S_temp:
      seqs_hash = copy.deepcopy(seqs_hash_temp)
      seqs_info = copy.deepcopy(seqs_info_temp)
      lowest_loss = loss
      lowest_loss_epoch = m+1

    else:
      loss = loss_i

"""## Grafico del Loss vs Epochs"""

x_axis = [x for x in range(epochs+1)]

fig = px.line(x=x_axis, y=loss_record, markers=True,
              labels={'x':'Epoch', 'y':'L(S)'})
fig.show()

"""## Set de primers finales"""

for index, seq_info in enumerate(seqs_info):
  print('Secuencia target:', seq_info['seq_id'])
  print('Forward = ',seq_info['forward'])
  print('Reverse = ',seq_info['reverse'])
  print('\n')

# import primer3
#
# final_primers =  []
#
# for index, seq_info in enumerate(seqs_info):
#   final_primers.append(primer3.calcTm(seq_info['forward']))
#   final_primers.append(primer3.calcTm(seq_info['reverse']))
#
# import plotly.express as px
# import numpy as np
#
# # Sample data
#
# fig = px.histogram(final_primers, nbins=50, title="Frequency Distribution Plot")
# fig.show()
#
# for index, seq_info in enumerate(seqs_info):
#   print('>primer_'+str(index)+'_f')
#   print(seq_info['forward'])
#   print('>primer_'+str(index)+'_r')
#   print(seq_info['reverse'])