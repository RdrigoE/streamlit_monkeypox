from numpy import info
import pandas as pd
import sys
import yaml


def get_info_from_structure(optional_df=None):
    # sample_name,pos,ref,alt,freq,type,mutation_anotation,apobec,depth
    if optional_df is None:
        df = pd.read_csv(sys.argv[1])
    else:
        df = optional_df
    info_dict = {}
# - Quantas variantes (validated) por amostra ?
# % number of variants per sample:
#     $ mean
#     $ median
#     $ std
    sample_variants = df.groupby('sample_name').size().describe()
    info_dict['total_entries'] = df.count()

    info_dict['total_entries']['total'] = info_dict['total_entries']['sample_name']
    info_dict['total_entries']['validated_variants'] = df[df['freq'] >= 0.51].count()[
        'sample_name']
    info_dict['total_entries']['minor_variants'] = df[df['freq'] < 0.51].count()[
        'sample_name']

    del info_dict['total_entries']['sample_name']
    del info_dict['total_entries']['pos']
    del info_dict['total_entries']['ref']
    del info_dict['total_entries']['alt']
    del info_dict['total_entries']['freq']
    del info_dict['total_entries']['type']
    del info_dict['total_entries']['depth']
    info_dict['all_variants'] = sample_variants
# All variants
    #print(sample_variants)

# % all validated variants
    validated_variants = df[df['freq'] >= 0.51].groupby(
        'sample_name').size().describe()
    #print(validated_variants)

    info_dict['validated_variants'] = validated_variants
# % all validated minor variants
    validated_minor_variants = df[df['freq'] < 0.51].groupby(
        'sample_name').size().describe()
    #print(validated_minor_variants)

    info_dict['validated_minor_variants'] = validated_minor_variants

# - Que tipo de variantes são (já fizeste isso, na discussão podes comentar
# se segue uma distribuição aleatória, e referir uma referência – a do APOBEC
# do Rambaud acho que refere isso)
#
# % number of types of variants

    total_variants = df.groupby('mutation_anotation').size()
    #print(total_variants)
    info_dict['total_mutations_annotation'] = total_variants
# % validated variants

    validated_variants_mut = df[df['freq'] >= 0.51].groupby(
        'mutation_anotation').size()
    #print(validated_variants_mut)

    info_dict['total_validated_mutations_annotation'] = validated_variants_mut
# % validated variants
    validated_minor_variants_mut = df[df['freq'] < 0.51].groupby(
        'mutation_anotation').size()
    #print(validated_minor_variants_mut)

    info_dict['validated_minor_variants_mut'] = validated_minor_variants_mut
# - relativamente à APOBEC, referir que é um padrão conhecido para alguns
# vírus como os Orthopox e mencionar a proporção que encontras, que está
# muito longe da distribuição aleatória [a referência do Rambaud novamente
# é boa para usar]
#
#
# % Find number of variants per sample with apobec pattern
#     $ mean
#     $ median
#     $ std

    apobec = df.groupby('apobec').size()
    #print(f"{apobec=}")

    info_dict['apobec'] = apobec
    apobec_val = df[df['freq'] >= 0.51].groupby('apobec').size()
    info_dict['apobec_validated'] = apobec_val
    #print(f"{apobec_val=}")

    apobec_minor = df[df['freq'] < 0.51].groupby('apobec').size()

    info_dict['apobec_minor'] = apobec_minor
    #print(f"{apobec_minor=}")
    for key in info_dict:
        info_dict[key] = info_dict[key].to_dict()
    return info_dict


if __name__ == "__main__":


    info_dict = get_info_from_structure()
    with open(sys.argv[2], "w") as handler:
        yaml.dump(info_dict, handler)
