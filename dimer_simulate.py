import primer3
import math
import time
import random
import primer3

 ############ Parametros termodinamicos deltaG NN ############

dg_5_overhangs_dict = { 'AA': -.51,
                        'AT': -.50,
                        'AC': -.96,
                        'AG': -.58,
                        'TA': -.71,
                        'TT': -.1,
                        'TC': -.58,
                        'TG': -.61,
                        'CA': -.42,
                        'CT': -.02,
                        'CC': -.52,
                        'CG': -.34,
                        'GA': -.62,
                        'GT': -.48,
                        'GC': -.72,
                        'GG': -.56,
                        }


dg_match_dict = {   'AA': -1,
                    'AT': -0.88,
                    'AC': -1.44,
                    'AG': -1.28,
                    'TA': -0.58,
                    'TT': -1,
                    'TC': -1.3,
                    'TG': -1.45,
                    'CA': -1.45,
                    'CT': -1.28,
                    'CC': -1.84,
                    'CG': -2.17,
                    'GA': -1.3,
                    'GT': -1.44,
                    'GC': -2.24,
                    'GG': -1.84,
                    }


dg_missmatchs_dict = {  'AG-T': 0.71,
                        'AT-G': 0.07,
                        'GT-G': -0.59,
                        'TG-T': 0.43,
                        'CA-C': 0.75,
                        'CC-A': 0.79,
                        'TC-A': 1.33,
                        'AC-T': 0.64,
                        'GC-T': 0.62,
                        'GT-C': 0.98,
                        'AG-A': 0.02,
                        'CA-G': 0.03,
                        'TA-G': 0.42,
                        'TG-A': 0.74,
                        'TA-A': 0.69,
                        'AC-C': 1.33,
                        'AG-G': -0.13,
                        'CG-G': -0.11,
                        'CT-T': -0.12,
                        'GT-T': 0.45,
                        'CG-T': -0.47,
                        'CT-G': -0.32,
                        'GG-T': 0.08,
                        'TT-G': 0.34,
                        'AA-C': 0.88,
                        'AC-A': 0.77,
                        'GA-C': 0.81,
                        'GC-A': 0.47,
                        'TA-C': 0.92,
                        'AT-C': 0.73,
                        'CC-T': 0.62,
                        'CT-C': 0.40,
                        'TC-T': 0.97,
                        'TT-C': 0.43, 
                        'AA-G': 0.14,
                        'CG-A': 0.11,
                        'GA-G': -0.25,
                        'GG-A': -0.52,
                        'AA-A': 0.61,
                        'CA-A': 0.43,
                        'GA-A': 0.17,
                        'CC-C': 0.70, 
                        'GC-C': 0.79, 
                        'TC-C': 1.05,
                        'GG-G': -1.11,
                        'TG-G': 0.44,
                        'AT-T': 0.69,
                        'TT-T': 0.68
                        } #48

                        

internal_loops_dict = { 3: 3.2,
                        4: 3.6,
                        5: 4.0,
                        6: 4.4,
                        7: 4.6,
                        8: 4.8,
                        }


double_mm_pen = 2




def complement_base(base):

    base = base.upper()

    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'

    else:
        raise KeyError

def score_alignment(long2mer, short2mer, double_mm_size):


    A2mer = long2mer
    B2mer = short2mer


    if len(A2mer) == 3 and len(B2mer) == 3:

        if (
            complement_base(A2mer[0]) != B2mer[0] and
            complement_base(A2mer[1]) != B2mer[1]
            ):

            double_mm_size += 1
            # print(A2mer[:2], '-', B2mer[:2])

            if (
                complement_base(A2mer[2]) == B2mer[2]
                ):

                # loop ends by match in 3rd base closing loop and is finally computed
                if double_mm_size  <= 6:
                    print('DoubleMM_size =', double_mm_size)
                    dg = internal_loops_dict[double_mm_size+2] #reemplazar con diccionario
                    double_mm_size = 0
                else:
                    dg = 4.6 + 2.44 * .002 * 310.15 * math.log((double_mm_size+2)/8)
                    print('DoubleMM_size =', double_mm_size)

                    double_mm_size = 0

            else:
                # loops continues in next round
                dg = 0


        elif (
            complement_base(A2mer[0]) == B2mer[0] and
            complement_base(A2mer[1]) != B2mer[1] and
            complement_base(A2mer[2]) == B2mer[2]
            ):
            # single mismatch

            mismatch_code_f = A2mer[:2] + '-' + B2mer[1]
            mismatch_code_r = B2mer[:0:-1] + '-' + A2mer[1]
            print('mismatch_code_f', mismatch_code_f)
            print('mismatch_code_r', mismatch_code_r)


            dg = dg_missmatchs_dict[mismatch_code_f]
            dg += dg_missmatchs_dict[mismatch_code_r]

            print('Second base mm')

        elif (
            complement_base(A2mer[0]) == B2mer[0] and
            complement_base(A2mer[1]) != B2mer[1] and
            complement_base(A2mer[2]) != B2mer[2]
            ):

            # start of double mismatch
            dg = 0

        elif (
            complement_base(A2mer[0]) == B2mer[0] and
            complement_base(A2mer[1]) == B2mer[1]
            ):
            dg = dg_match_dict[A2mer[:2]]
            print('Match')



        elif (complement_base(A2mer[0]) != B2mer[0] and
            complement_base(A2mer[1]) == B2mer[1]
            ):

            # only first base mismatch, ignore.
            print('First base mm (ignore)')

            dg = 0

        else:
            raise KeyError

    elif len(A2mer) == 2 or len(B2mer) == 2:


        if (
            complement_base(A2mer[0]) != B2mer[0] and
            complement_base(A2mer[1]) != B2mer[1]
            ):

            double_mm_size += 1

            # loop ends by end of oligo and is finally computed

            if double_mm_size  <= 6:
                print('size =', double_mm_size)
                dg = internal_loops_dict[double_mm_size+2]
                double_mm_size = 0

            else:
                dg = 4.6 + 2.44 * .002 * 310.15 * math.log((double_mm_size+2)/8)

                double_mm_size = 0



        elif (
            complement_base(A2mer[0]) == B2mer[0] and
            complement_base(A2mer[1]) != B2mer[1]
            ):
            # single mismatch (2nd base)

            mismatch_code = A2mer[:2] + '-' + B2mer[1]

            dg = dg_missmatchs_dict[mismatch_code]


        elif (
            complement_base(A2mer[0]) == B2mer[0] and
            complement_base(A2mer[1]) == B2mer[1]
            ):
            dg = dg_match_dict[A2mer[:2]]

            if len(A2mer) == 2 and len(B2mer) == 3:
                # Overhang
                dg += dg_5_overhangs_dict[B2mer[::-1][:2]]



        elif (complement_base(A2mer[0]) != B2mer[0] and
            complement_base(A2mer[1]) == B2mer[1]
            ):
            # match in end ignored???????

            if len(A2mer) == 2 and len(B2mer) == 3:
                # Overhang
                dg = dg_5_overhangs_dict[B2mer[::-1][:2]]

            else:
                # first base mismatch, ignore.
                dg = 0

        else:
            raise KeyError        

    return dg, double_mm_size


def dimer_simulate(primerA, primerB):

    len_A = len(primerA)
    len_B = len(primerB)

    min_dg = 0

    if len_A > len_B:
        long_primer = primerA
        short_primer = primerB[::-1]

    else:
        long_primer = primerB
        short_primer = primerA[::-1]

    n_2mers = len(short_primer) - 1

    dg = 0

    for pos_index in range(len(long_primer)):
        print('ESTRUCTURA', str(pos_index+1))
        print(long_primer)
        print(' '*pos_index+short_primer)
        dg = 0
        double_mm_size = 0

        for short_index in range(n_2mers):
            short_2mer = short_primer[short_index:short_index+3]
            long_index = pos_index+short_index
            long_2mer = long_primer[long_index:long_index+3]
            print(short_2mer[:2], '-', long_2mer[:2])


            if len(short_2mer) > 1 and len(long_2mer) > 1:
                dg_2mer, double_mm_size = score_alignment(long_2mer, short_2mer, double_mm_size)
                dg += dg_2mer
                print('ΔG 2mer =', dg_2mer)
            else:
                break

            # print('ΔG total =', dg)
            if dg < min_dg:
                min_dg = round(dg, 3)
                min_pos = pos_index

        print('ΔG_final' +' = '+ str(dg))


    return min_dg, min_pos


def random_seq(size):
    seq = ''
    for _ in range(size):
        base = random.choice(['A', 'C', 'T', 'G'])
        seq += base

    return ''.join(seq)


# seqs=5000
# min_dg = 50
start = time.time()
# for i in range(seqs):
#     random_primer1 = random_seq(25)
#     random_primer2 = random_seq(23)
#     dg, pos = dimer_simulate(random_primer1, random_primer2)
#     if dg < min_dg:
#         best_1 = random_primer1
#         best_2 = random_primer2
#         min_dg = dg
#         min_pos = pos

end = time.time()




print('tiempo =', end-start)
# print('tiempo/seq =', (end-start)/seqs)


min_dg = 100

primer1 = 'GGACTGACG'.replace(' ', '')
primer2 = 'CGTCGGTCC'.replace(' ', '')

dg, pos = dimer_simulate(primer1, primer2)
if dg < min_dg:
    best_1 = primer1
    best_2 = primer2
    min_dg = dg
    min_pos = pos


print('\n\n\n\nDimero mas estable\n')
if len(primer1) > len(primer2):
    print(primer1)
    print(' '*min_pos+primer2[::-1])
    print('\n\nprimer1 = ', primer1)
    print('primer2 = ', primer2)

else:
    print(primer2)
    print(' '*min_pos+primer1[::-1])
    print('\n\nprimer1 = ', primer2)
    print('primer2 = ', primer1)

print('ΔG =', min_dg)


print(primer3.calcHeterodimer(primer1, primer2))
dimer = primer3.calcHeterodimer(primer1, primer2, output_structure=True)

for i in dimer.ascii_structure_lines:
    print(i)