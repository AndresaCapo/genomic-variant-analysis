import argparse
import csv

import pysam

from analise_variantes import GENES, obter_cobertura
from suporte_alelico import contar_snp, contar_delecao


COLUNAS = [
    "gene", "chrom", "pos", "ref", "alt",
    "tipo", "gt", "filter", "dp_vcf", "ad_vcf",
    "dp_bam_inicial", "status",
    "ref_alinhamentos", "alt_alinhamentos",
    "outro_alinhamentos", "nao_classificado_alinhamentos",
    "ref_fragmentos", "alt_fragmentos",
    "discordante_fragmentos", "nao_classificado_fragmentos",
]


def classificar_tipo(ref, alts):
    if not alts:
        return "SEM_ALT"

    if len(alts) > 1:
        return "MULTIALELICA"

    alt = alts[0]

    if alt == "*" or alt.startswith("<") or "[" in alt or "]" in alt:
        return "ALELO_SIMBOLICO"

    if len(ref) == len(alt) == 1:
        return "SNP"

    if len(ref) == len(alt):
        return "MNV"

    if len(ref) > 1 and len(alt) == 1 and ref.startswith(alt):
        return "DELECAO_SIMPLES"

    if len(ref) == 1 and len(alt) > 1 and alt.startswith(ref):
        return "INSERCAO_SIMPLES"

    return "VARIANTE_COMPLEXA"


def formatar_valores(valores):
    if valores is None:
        return ""
    return ",".join("." if x is None else str(x) for x in valores)


def exportar(vcf_path, bam_path, saida):
    total = 0

    with pysam.VariantFile(vcf_path) as vcf, \
         pysam.AlignmentFile(bam_path, "rb") as bam, \
         open(saida, "w", newline="", encoding="utf-8") as arquivo:

        escritor = csv.DictWriter(arquivo, fieldnames=COLUNAS)
        escritor.writeheader()

        for gene, (chrom, inicio, fim) in GENES.items():
            for variante in vcf.fetch(chrom, inicio - 1, fim):
                if not inicio <= variante.pos <= fim:
                    continue

                amostra = next(iter(variante.samples.values()), None)
                gt = amostra.get("GT") if amostra is not None else None

                if gt is None:
                    genotipo = "."
                else:
                    sep = "|" if amostra.phased else "/"
                    genotipo = sep.join(
                        "." if x is None else str(x) for x in gt
                    )

                alts = variante.alts or ()
                dp_bam = obter_cobertura(bam, variante)

                linha = {
                    "gene": gene,
                    "chrom": variante.contig,
                    "pos": variante.pos,
                    "ref": variante.ref,
                    "alt": ",".join(alts),
                    "tipo": classificar_tipo(variante.ref, alts),
                    "gt": genotipo,
                    "filter": ",".join(variante.filter.keys()) or ".",
                    "dp_vcf": amostra.get("DP") if amostra is not None else "",
                    "ad_vcf": formatar_valores(
                        amostra.get("AD") if amostra is not None else None
                    ),
                    "dp_bam_inicial": dp_bam,
                }

                if dp_bam == 0:
                    linha["status"] = "SEM_COBERTURA_INICIAL"

                elif len(alts) != 1:
                    linha["status"] = "TIPO_NAO_SUPORTADO"

                else:
                    alt = alts[0]

                    if len(variante.ref) == 1 and len(alt) == 1:
                        linha["tipo"] = "SNP"
                        resultado = contar_snp(
                            bam, variante.contig, variante.pos,
                            variante.ref, alt
                        )

                    elif (
                        len(variante.ref) > 1
                        and len(alt) == 1
                        and variante.ref.startswith(alt)
                    ):
                        linha["tipo"] = "DELECAO_SIMPLES"
                        resultado = contar_delecao(
                            bam, variante.contig, variante.pos,
                            variante.ref, alt
                        )

                    else:
                        resultado = None

                    if resultado is None:
                        linha["status"] = "TIPO_NAO_SUPORTADO"
                    else:
                        linha["status"] = "SUPORTE_AVALIADO"

                        a = resultado["alinhamentos"]
                        f = resultado["fragmentos"]

                        linha.update({
                            "ref_alinhamentos": a.get("REF", 0),
                            "alt_alinhamentos": a.get("ALT", 0),
                            "ref_fragmentos": f.get("REF", 0),
                            "alt_fragmentos": f.get("ALT", 0),
                            "discordante_fragmentos": f.get("DISCORDANTE", 0),
                        })

                        if linha["tipo"] == "SNP":
                            linha["outro_alinhamentos"] = a.get("OUTRO", 0)
                        else:
                            linha["nao_classificado_alinhamentos"] = (
                                a.get("NAO_CLASSIFICADO", 0)
                            )
                            linha["nao_classificado_fragmentos"] = (
                                f.get("NAO_CLASSIFICADO", 0)
                            )

                escritor.writerow(linha)
                total += 1

    print(f"Arquivo gerado: {saida}")
    print(f"Total de variantes exportadas: {total}")


def main():
    parser = argparse.ArgumentParser(
        description="Exporta resultados da analise genetica para CSV"
    )
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--saida", required=True)

    args = parser.parse_args()
    exportar(args.vcf, args.bam, args.saida)


if __name__ == "__main__":
    main()
