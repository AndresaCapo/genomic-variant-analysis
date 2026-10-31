import argparse

import pysam
from suporte_alelico import contar_snp, contar_delecao


GENES = {
    "TGFBR2": ("chr3", 30606356, 30694142),
    "HEXA": ("chr15", 72340924, 72376014),
    "FANCA": ("chr16", 89737549, 89816647),
}


def obter_cobertura(bam, variante):
    cobertura = bam.count_coverage(
        variante.contig,
        variante.pos - 1,
        variante.pos,
        quality_threshold=0,
        read_callback="all",
    )
    return sum(bases[0] for bases in cobertura)



def analisar(vcf_path, bam_path):
    with pysam.VariantFile(vcf_path) as vcf, \
         pysam.AlignmentFile(bam_path, "rb") as bam:

        for gene, (chrom, inicio, fim) in GENES.items():
            total = 0
            cobertas = 0
            nao_suportadas = 0

            print(f"\nGene: {gene}")

            for variante in vcf.fetch(chrom, inicio - 1, fim):
                if not inicio <= variante.pos <= fim:
                    continue

                total += 1
                dp = obter_cobertura(bam, variante)

                alts = variante.alts or ()

                descricao = (
                    f"{variante.contig}:{variante.pos} "
                    f"{variante.ref}>{','.join(alts) if alts else '.'}"
                )

                amostras = list(variante.samples.values())

                if amostras:
                    gt = amostras[0].get("GT")

                    if gt is not None:
                        separador = "|" if amostras[0].phased else "/"
                        genotipo = separador.join(
                            "." if alelo is None else str(alelo)
                            for alelo in gt
                        )
                    else:
                        genotipo = "."
                else:
                    genotipo = "."

                filtros = ",".join(variante.filter.keys()) or "."

                print(f"  Variante: {descricao}")
                print(f"  GT: {genotipo} | FILTER: {filtros}")
                print(f"  Cobertura na posicao inicial: {dp}")

                if dp == 0:
                    print("  Status: SEM_COBERTURA_NA_POSICAO_INICIAL")
                    continue

                cobertas += 1

                if len(alts) != 1:
                    resultado = None
                    tipo = "Multialelica ou sem ALT"

                elif len(variante.ref) == 1 and len(alts[0]) == 1:
                    tipo = "SNP"
                    resultado = contar_snp(
                        bam, variante.contig, variante.pos,
                        variante.ref, alts[0]
                    )

                elif (
                    len(variante.ref) > 1
                    and len(alts[0]) == 1
                    and variante.ref.startswith(alts[0])
                ):
                    tipo = "Delecao simples"
                    resultado = contar_delecao(
                        bam, variante.contig, variante.pos,
                        variante.ref, alts[0]
                    )

                else:
                    tipo = "Tipo ainda nao suportado"
                    resultado = None

                print(f"  Tipo: {tipo}")

                if resultado is None:
                    nao_suportadas += 1
                    print("  Suporte alelico: NAO AVALIADO")
                    continue

                a = resultado["alinhamentos"]
                f = resultado["fragmentos"]

                print(
                    "  Alinhamentos: "
                    f"REF={a.get('REF', 0)} "
                    f"ALT={a.get('ALT', 0)} "
                    f"OUTRO={a.get('OUTRO', 0)} "
                    f"NAO_CLASSIFICADO={a.get('NAO_CLASSIFICADO', 0)}"
                )

                print(
                    "  Fragmentos: "
                    f"REF={f.get('REF', 0)} "
                    f"ALT={f.get('ALT', 0)} "
                    f"DISCORDANTE={f.get('DISCORDANTE', 0)} "
                    f"NAO_CLASSIFICADO={f.get('NAO_CLASSIFICADO', 0)}"
                )

            print(f"  Total de variantes no VCF: {total}")
            print(f"  Com cobertura inicial: {cobertas}")
            print(f"  Sem cobertura inicial: {total - cobertas}")
            print(f"  Tipos ainda nao suportados: {nao_suportadas}")

def main():
    parser = argparse.ArgumentParser(
        description="Analise de variantes em TGFBR2, HEXA e FANCA (GRCh38)"
    )
    parser.add_argument("--vcf", required=True, help="VCF compactado e indexado")
    parser.add_argument("--bam", required=True, help="BAM indexado")

    args = parser.parse_args()
    analisar(args.vcf, args.bam)


if __name__ == "__main__":
    main()
