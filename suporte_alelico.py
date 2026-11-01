"""Avaliacao de suporte alelico em alinhamentos de sequenciamento."""

from collections import Counter, defaultdict


def contar_snp(bam, chrom, pos, ref, alt, min_mapq=20, min_baseq=20):
    """Conta o suporte REF e ALT para um SNP bialelico.
    Aplica filtros MAPQ e BaseQ e apresenta contagens por alinhamento
    e por QNAME. Observacoes excluidas podem nao integrar as contagens."""
    if len(ref) != 1 or len(alt) != 1:
        raise ValueError("Esta funcao aceita apenas SNPs bialelicos.")

    alinhamentos = Counter()
    fragmentos = defaultdict(set)

    for read in bam.fetch(chrom, pos - 1, pos):
        if (
            read.is_unmapped
            or read.is_duplicate
            or read.is_secondary
            or read.is_supplementary
            or read.is_qcfail
            or read.mapping_quality < min_mapq
        ):
            continue

        if read.query_sequence is None or read.query_qualities is None:
            continue

        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            if rpos != pos - 1:
                continue

            if read.query_qualities[qpos] < min_baseq:
                break

            base = read.query_sequence[qpos].upper()

            if base == ref.upper():
                alelo = "REF"
            elif base == alt.upper():
                alelo = "ALT"
            else:
                alelo = "OUTRO"

            alinhamentos[alelo] += 1
            fragmentos[read.query_name].add(alelo)
            break

    contagem_fragmentos = Counter()

    for evidencias in fragmentos.values():
        if len(evidencias) == 1:
            contagem_fragmentos[next(iter(evidencias))] += 1
        else:
            contagem_fragmentos["DISCORDANTE"] += 1

    return {
        "alinhamentos": dict(alinhamentos),
        "fragmentos": dict(contagem_fragmentos),
    }


def contar_delecao(bam, chrom, pos, ref, alt, min_mapq=20, min_baseq=20):
    """Avalia o suporte a uma delecao simples bialelica.
    Examina a sequencia e as operacoes CIGAR dos alinhamentos,
    considerando filtros de qualidade e a base ancora.
    Exige delecao no CIGAR na posicao e no tamanho esperados.
    Distingue suporte REF, ALT e evidencias nao classificaveis."""
    if len(alt) != 1 or len(ref) <= 1 or not ref.startswith(alt):
        raise ValueError("A funcao aceita apenas delecoes simples bialelicas.")

    inicio = pos - 1
    fim = inicio + len(ref)
    tamanho_delecao = len(ref) - 1

    alinhamentos = Counter()
    fragmentos = defaultdict(set)

    for read in bam.fetch(chrom, inicio, fim):
        if (
            read.is_unmapped
            or read.is_duplicate
            or read.is_secondary
            or read.is_supplementary
            or read.is_qcfail
            or read.mapping_quality < min_mapq
        ):
            continue

        if read.query_sequence is None or read.query_qualities is None:
            continue

        # Preserva QNAMEs sem evidencia classificavel.
        fragmentos[read.query_name]

        mapa = {
            r: q
            for q, r in read.get_aligned_pairs(matches_only=False)
            if r is not None and inicio <= r < fim
        }

        def base_valida(coordenada):
            q = mapa.get(coordenada)
            if q is None or read.query_qualities[q] < min_baseq:
                return None
            return read.query_sequence[q].upper()

        if base_valida(inicio) != alt.upper():
            alinhamentos["NAO_CLASSIFICADO"] += 1
            continue

        referencia = all(
            base_valida(inicio + i) == base.upper()
            for i, base in enumerate(ref)
        )

        cursor = read.reference_start
        delecao = False

        for op, tamanho in read.cigartuples or []:
            if (
                op == 2
                and cursor == inicio + 1
                and tamanho == tamanho_delecao
            ):
                delecao = True

            # Apenas operacoes que consomem referencia avancam o cursor.
            if op in (0, 2, 3, 7, 8):
                cursor += tamanho

        if referencia:
            alelo = "REF"
        elif delecao:
            alelo = "ALT"
        else:
            alinhamentos["NAO_CLASSIFICADO"] += 1
            continue

        alinhamentos[alelo] += 1
        fragmentos[read.query_name].add(alelo)

    contagem_fragmentos = Counter()

    for evidencias in fragmentos.values():
        if not evidencias:
            contagem_fragmentos["NAO_CLASSIFICADO"] += 1
        elif len(evidencias) == 1:
            contagem_fragmentos[next(iter(evidencias))] += 1
        else:
            contagem_fragmentos["DISCORDANTE"] += 1

    return {
        "alinhamentos": dict(alinhamentos),
        "fragmentos": dict(contagem_fragmentos),
    }
