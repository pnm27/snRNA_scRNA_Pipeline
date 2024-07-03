# Mem allocation per thread (here thread is core/tasks), everything is in MB
def allocate_mem_SS(wildcards, attempt):
    return 25000+1000*(attempt-1)


def allocate_mem_PGCB(wildcards, attempt):
    return 2500*attempt+2500


def allocate_mem_PRNA(wildcards, attempt):
    return 2500*attempt+2500


def allocate_mem_CFB(wildcards, attempt):
    return 50*attempt+80


def allocate_mem_CMF(wildcards, attempt):
    return 50*attempt+70


def allocate_mem_BKI(wildcards, attempt):
    return 1000*attempt+1000


def allocate_mem_RK(wildcards, attempt):
    return 1000*attempt + 500


def allocate_mem_RBCor(wildcards, attempt):
    return 1500*attempt+1500


def allocate_mem_RBS(wildcards, attempt):
    return 3000*attempt+3000


def allocate_mem_RBCnt(wildcards, attempt):
    return 1500*attempt+1500


def allocate_mem_CHB(wildcards, attempt):
    return 1500*attempt+1500


def allocate_mem_RCS(wildcards, attempt):
    return 4500*attempt+2500


def allocate_mem_DXP(wildcards, attempt):
    return 3500*attempt+3500


def allocate_mem_CICS(wildcards, attempt):
    return 2000+500*(attempt-1)


def allocate_mem_CISPD(wildcards, attempt):
    return 500+500*(attempt-1)


def allocate_mem_GIH(wildcards, attempt):
    return 3500*attempt+3500 # Not finalised yet


def allocate_mem_cS(wildcards, attempt):
    if 'vcf_type' in wildcards:
        if wildcards.vcf_type.endswith('_max50k'):
            return 200+80*(attempt-1)
        elif wildcards.vcf_type.endswith('25pct'):
            return 800+100*(attempt-1)
        elif wildcards.vcf_type.endswith('gencode'):
            return 1300+200*(attempt-1)
        else:
            return 1500+200*(attempt-1)
    else:
        return 1500+200*(attempt-1)


def allocate_mem_vS(wildcards, attempt):
    if 'vcf_type' in wildcards:
        if wildcards.vcf_type.endswith('_max50k'):
            return 3000+100*(attempt-1)
        elif wildcards.vcf_type.endswith('25pct'):
            return 5000+400*(attempt-1)
        elif wildcards.vcf_type.endswith('gencode'):
            return 10000+600*(attempt-1)
        else:
            return 22000+500*(attempt-1)
    else:
        return 22000+500*(attempt-1)


def allocate_mem_AOTFCM(wildcards, attempt):
    return 2000+500*(attempt-1)

def allocate_mem_US(wildcards, attempt):
    return 250*attempt+250


def allocate_mem_CIS(wildcards, attempt):
    return 500*attempt+500


def allocate_mem_FCB(wildcards, attempt):
    return 150*attempt+150


def allocate_mem_SB(wildcards, attempt):
    return 75*attempt+75


def allocate_mem_SSSA(wildcards, attempt):
    return 6000+1000*(attempt-1)


def allocate_mem_CPDV(wildcards, attempt):
    return 75*attempt+75


def allocate_mem_SDB(wildcards, attempt):
    return 2500


def allocate_mem_RN(wildcards, attempt):
    return 20000+10000*(attempt-1)


def allocate_mem_QM(wildcards):
    return 5000


def allocate_mem_MBS(wildcards, attempt):
    return 400 + 200*(attempt-1)


def allocate_mem_PC(wildcards, attempt):
    return 200 + 500*(attempt-1)


# Time allocation in minutes (doesn't include for the rules cellSNP and vireoSNP)
# CAREFULL about exceeding time limits
def allocate_time_SS(wildcards, attempt):
    return 1440


def allocate_time_PGCB(wildcards, attempt):
    return 60*attempt+60


def allocate_time_PRNA(wildcards, attempt):
    return 210*attempt+210


def allocate_time_CFB(wildcards, attempt):
    return 5*attempt+5


def allocate_time_CMF(wildcards, attempt):
    return 5*attempt+5


def allocate_time_BKI(wildcards, attempt):
    return 5*attempt+5


def allocate_time_RK(wildcards, attempt):
    return 20*attempt + 20


def allocate_time_RBCor(wildcards, attempt):
    return 5*attempt+5


def allocate_time_RBS(wildcards, attempt):
    return 5*attempt+5


def allocate_time_RBCnt(wildcards, attempt):
    return 5*attempt+5


def allocate_time_CHB(wildcards, attempt):
    return 5*attempt+5


def allocate_time_RCS(wildcards, attempt):
    return 5*attempt+5


def allocate_time_DXP(wildcards, attempt):
    return 15*attempt+15


def allocate_time_CICS(wildcards, attempt):
    return 2*attempt+1


def allocate_time_GIH(wildcards, attempt):
    return 3500*attempt+3500 # Not finalised yet


def allocate_time_cS(wildcards, attempt):
    if 'vcf_type' in wildcards:
        if wildcards.vcf_type.endswith('_max50k'):
            return 30+20*(attempt-1)
        elif wildcards.vcf_type.endswith('25pct'):
            return 70+20*(attempt-1)
        elif wildcards.vcf_type.endswith('gencode'):
            return 180+40*(attempt-1)
        else:
            return 240+60*(attempt-1)
    else:
        return 240+60*(attempt-1)


def allocate_time_vS(wildcards, attempt):
    if 'vcf_type' in wildcards:
        if wildcards.vcf_type.endswith('_max50k'):
            return 10+5*(attempt-1)
        elif wildcards.vcf_type.endswith('25pct'):
            return 70+20*(attempt-1)
        elif wildcards.vcf_type.endswith('gencode'):
            return 180+40*(attempt-1)
        else:
            return 180+60*(attempt-1)
    else:
        return 180+60*(attempt-1)


def allocate_time_AOTFCM(wildcards, attempt):
    return 15+10*(attempt-1)


def allocate_time_US(wildcards, attempt):
    return 15*attempt+15


def allocate_time_FCB(wildcards, attempt):
    return 20*attempt+20


def allocate_time_SSSA(wildcards, attempt):
    return 120


def allocate_time_CPDV(wildcards, attempt):
    return 300


def allocate_time_SDB(wildcards, attempt):
    return 180


def allocate_timem_RN(wildcards, attempt):
    return 360


def allocate_time_QM(wildcards, attempt):
    return 60 + 20*(attempt-1)


def allocate_time_MBS(wildcards, attempt):
    return 300 + 30*(attempt-1)


def allocate_time_PC(wildcards, attempt):
    return 10 + 200*(attempt-1)

