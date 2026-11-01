"""Testes automatizados com dados sinteticos para validar o projeto."""

import tempfile
import unittest
from pathlib import Path

import pysam

from exportar_resultados import classificar_tipo, avaliar_confirmacao_bam, COLUNAS
from suporte_alelico import contar_snp, contar_delecao


def criar_bam(leituras):
    """Cria um BAM sintetico indexado para testes de suporte alelico."""
    pasta = tempfile.TemporaryDirectory()
    caminho = str(Path(pasta.name) / "teste.bam")

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrTeste", "LN": 1000}],
    }

    with pysam.AlignmentFile(caminho, "wb", header=header) as bam:
        for nome, sequencia, cigar, qualidade in leituras:
            read = pysam.AlignedSegment()
            read.query_name = nome
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 99
            read.mapping_quality = 60
            read.query_sequence = sequencia
            read.cigartuples = cigar
            read.query_qualities = pysam.qualitystring_to_array(qualidade)
            bam.write(read)

    pysam.index(caminho)
    return pasta, caminho


class TestesRegressao(unittest.TestCase):

    def test_classificacao_snp(self):
        self.assertEqual(classificar_tipo("A", ("G",)), "SNP")

    def test_classificacao_delecao(self):
        self.assertEqual(
            classificar_tipo("CAG", ("C",)),
            "DELECAO_SIMPLES"
        )

    def test_classificacao_insercao(self):
        self.assertEqual(
            classificar_tipo("C", ("CT",)),
            "INSERCAO_SIMPLES"
        )

    def test_contagem_snp(self):
        leituras = [
            ("frag1", "A", [(0, 1)], "I"),
            ("frag1", "A", [(0, 1)], "I"),
            ("frag2", "G", [(0, 1)], "I"),
        ]

        pasta, caminho = criar_bam(leituras)

        try:
            with pysam.AlignmentFile(caminho, "rb") as bam:
                r = contar_snp(bam, "chrTeste", 100, "A", "G")

            self.assertEqual(r["alinhamentos"]["REF"], 2)
            self.assertEqual(r["alinhamentos"]["ALT"], 1)
            self.assertEqual(r["fragmentos"]["REF"], 1)
            self.assertEqual(r["fragmentos"]["ALT"], 1)
        finally:
            pasta.cleanup()

    def test_contagem_delecao(self):
        leituras = [
            ("frag1", "CAG", [(0, 3)], "III"),
            ("frag2", "C", [(0, 1), (2, 2)], "I"),
            ("frag3", "C", [(0, 1)], "I"),
        ]

        pasta, caminho = criar_bam(leituras)

        try:
            with pysam.AlignmentFile(caminho, "rb") as bam:
                r = contar_delecao(
                    bam, "chrTeste", 100, "CAG", "C"
                )

            self.assertEqual(r["alinhamentos"]["REF"], 1)
            self.assertEqual(r["alinhamentos"]["ALT"], 1)
            self.assertEqual(
                r["alinhamentos"]["NAO_CLASSIFICADO"], 1
            )
        finally:
            pasta.cleanup()



class TestesAdicionais(unittest.TestCase):

    def test_snp_baixa_qualidade(self):
        leituras = [
            ("frag1", "A", [(0, 1)], "I"),
            ("frag2", "G", [(0, 1)], "!"),
        ]
        pasta, caminho = criar_bam(leituras)
        try:
            with pysam.AlignmentFile(caminho, "rb") as bam:
                r = contar_snp(bam, "chrTeste", 100, "A", "G")
            self.assertEqual(r["alinhamentos"].get("REF", 0), 1)
            self.assertEqual(r["alinhamentos"].get("ALT", 0), 0)
        finally:
            pasta.cleanup()

    def test_snp_fragmento_discordante(self):
        leituras = [
            ("frag1", "A", [(0, 1)], "I"),
            ("frag1", "G", [(0, 1)], "I"),
        ]
        pasta, caminho = criar_bam(leituras)
        try:
            with pysam.AlignmentFile(caminho, "rb") as bam:
                r = contar_snp(bam, "chrTeste", 100, "A", "G")
            self.assertEqual(r["fragmentos"]["DISCORDANTE"], 1)
        finally:
            pasta.cleanup()

    def test_delecao_com_softclip(self):
        leituras = [
            ("frag1", "TTC", [(4, 2), (0, 1), (2, 2)], "III"),
        ]
        pasta, caminho = criar_bam(leituras)
        try:
            with pysam.AlignmentFile(caminho, "rb") as bam:
                r = contar_delecao(bam, "chrTeste", 100, "CAG", "C")
            self.assertEqual(r["alinhamentos"]["ALT"], 1)
        finally:
            pasta.cleanup()

    def test_delecao_baixa_qualidade(self):
        leituras = [
            ("frag1", "C", [(0, 1), (2, 2)], "!"),
        ]
        pasta, caminho = criar_bam(leituras)
        try:
            with pysam.AlignmentFile(caminho, "rb") as bam:
                r = contar_delecao(bam, "chrTeste", 100, "CAG", "C")
            self.assertEqual(
                r["alinhamentos"]["NAO_CLASSIFICADO"], 1
            )
        finally:
            pasta.cleanup()

    def test_delecao_fragmento_discordante(self):
        leituras = [
            ("frag1", "CAG", [(0, 3)], "III"),
            ("frag1", "C", [(0, 1), (2, 2)], "I"),
        ]
        pasta, caminho = criar_bam(leituras)
        try:
            with pysam.AlignmentFile(caminho, "rb") as bam:
                r = contar_delecao(bam, "chrTeste", 100, "CAG", "C")
            self.assertEqual(r["fragmentos"]["DISCORDANTE"], 1)
        finally:
            pasta.cleanup()



class TestesConfirmacaoBam(unittest.TestCase):

    def avaliar(self, **alteracoes):
        linha = {
            "status": "SUPORTE_AVALIADO",
            "filter": "PASS",
            "gt": "0/1",
            "ref_fragmentos": 0,
            "alt_fragmentos": 0,
        }
        linha.update(alteracoes)
        return avaliar_confirmacao_bam(linha)

    def test_heterozigoto_compativel(self):
        self.assertEqual(len(COLUNAS), 21)
        self.assertEqual(COLUNAS[-1], "confirmacao_bam")
        self.assertEqual(
            self.avaliar(ref_fragmentos=20, alt_fragmentos=18),
            "SUPORTE_COMPATIVEL_COM_GENOTIPO"
        )

    def test_heterozigoto_faseado(self):
        self.assertEqual(
            self.avaliar(gt="0|1", ref_fragmentos=10, alt_fragmentos=12),
            "SUPORTE_COMPATIVEL_COM_GENOTIPO"
        )

    def test_heterozigoto_limite_minimo(self):
        self.assertEqual(
            self.avaliar(ref_fragmentos=3, alt_fragmentos=3),
            "SUPORTE_COMPATIVEL_COM_GENOTIPO"
        )

    def test_heterozigoto_suporte_insuficiente(self):
        self.assertEqual(
            self.avaliar(ref_fragmentos=20, alt_fragmentos=2),
            "SUPORTE_INSUFICIENTE"
        )

    def test_homozigoto_alternativo_compativel(self):
        self.assertEqual(
            self.avaliar(gt="1/1", ref_fragmentos=0, alt_fragmentos=25),
            "SUPORTE_COMPATIVEL_COM_GENOTIPO"
        )

    def test_homozigoto_alternativo_abaixo_90_porcento(self):
        self.assertEqual(
            self.avaliar(gt="1/1", ref_fragmentos=12, alt_fragmentos=100),
            "SUPORTE_INSUFICIENTE"
        )

    def test_sem_cobertura(self):
        self.assertEqual(
            self.avaliar(status="SEM_COBERTURA_INICIAL"),
            "NAO_AVALIAVEL_SEM_COBERTURA"
        )

    def test_tipo_nao_suportado(self):
        self.assertEqual(
            self.avaliar(status="TIPO_NAO_SUPORTADO"),
            "NAO_AVALIAVEL_TIPO"
        )

    def test_filtro_vcf_nao_pass(self):
        self.assertEqual(
            self.avaliar(
                filter="FBQualDepth",
                ref_fragmentos=15,
                alt_fragmentos=15
            ),
            "NAO_AVALIAVEL_FILTRO_VCF"
        )

    def test_genotipo_fora_do_escopo(self):
        self.assertEqual(
            self.avaliar(gt="0/0", ref_fragmentos=20),
            "NAO_AVALIAVEL_GENOTIPO"
        )



class TesteIntegracao(unittest.TestCase):

    def test_exportacao_completa_vcf_bam_csv(self):
        import csv
        from exportar_resultados import exportar

        pos_snp = 30671634
        pos_sem_cobertura = 30671650
        pos_insercao = 30671670

        with tempfile.TemporaryDirectory() as pasta:
            diretorio = Path(pasta)
            vcf_txt = diretorio / "sintetico.vcf"
            vcf_gz = diretorio / "sintetico.vcf.gz"
            bam_path = diretorio / "sintetico.bam"
            csv_path = diretorio / "resultado.csv"

            header = {
                "HD": {"VN": "1.6", "SO": "coordinate"},
                "SQ": [
                    {"SN": "chr3", "LN": 100000000},
                    {"SN": "chr15", "LN": 100000000},
                    {"SN": "chr16", "LN": 100000000},
                ],
            }

            with pysam.AlignmentFile(
                str(bam_path), "wb", header=header
            ) as bam:
                leituras = (
                    [(pos_snp, "A")] * 3
                    + [(pos_snp, "G")] * 3
                    + [(pos_insercao, "C")] * 2
                )

                for indice, (posicao, base) in enumerate(leituras):
                    leitura = pysam.AlignedSegment()
                    leitura.query_name = f"frag_{indice}"
                    leitura.flag = 0
                    leitura.reference_id = 0
                    leitura.reference_start = posicao - 1
                    leitura.mapping_quality = 60
                    leitura.cigartuples = [(0, 1)]
                    leitura.query_sequence = base
                    leitura.query_qualities = (
                        pysam.qualitystring_to_array("I")
                    )
                    bam.write(leitura)

            pysam.index(str(bam_path))

            cabecalho_vcf = [
                "##fileformat=VCFv4.2",
                "##contig=<ID=chr3,length=100000000>",
                "##contig=<ID=chr15,length=100000000>",
                "##contig=<ID=chr16,length=100000000>",
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tAMOSTRA",
            ]

            registros = [
                f"chr3\t{pos_snp}\t.\tA\tG\t.\tPASS\t.\tGT:DP:AD\t0/1:6:3,3",
                f"chr3\t{pos_sem_cobertura}\t.\tC\tT\t.\tPASS\t.\tGT:DP:AD\t0/1:0:0,0",
                f"chr3\t{pos_insercao}\t.\tC\tCA\t.\tPASS\t.\tGT:DP:AD\t0/1:2:1,1",
                "chr15\t1000\t.\tA\tG\t.\tPASS\t.\tGT:DP:AD\t0/1:0:0,0",
                "chr16\t1000\t.\tA\tG\t.\tPASS\t.\tGT:DP:AD\t0/1:0:0,0",
            ]

            vcf_txt.write_text(
                "\n".join(cabecalho_vcf + registros) + "\n",
                encoding="utf-8"
            )

            pysam.tabix_compress(
                str(vcf_txt), str(vcf_gz), force=True
            )
            pysam.tabix_index(
                str(vcf_gz), preset="vcf", force=True
            )

            exportar(
                str(vcf_gz),
                str(bam_path),
                str(csv_path)
            )

            with csv_path.open(
                newline="", encoding="utf-8"
            ) as resultado:
                leitor = csv.DictReader(resultado)
                colunas = leitor.fieldnames
                linhas = list(leitor)

            self.assertEqual(colunas, COLUNAS)
            self.assertEqual(len(colunas), 21)
            self.assertEqual(len(linhas), 3)

            por_posicao = {
                int(linha["pos"]): linha
                for linha in linhas
            }

            snp = por_posicao[pos_snp]
            self.assertEqual(snp["gene"], "TGFBR2")
            self.assertEqual(snp["status"], "SUPORTE_AVALIADO")
            self.assertEqual(snp["ref_fragmentos"], "3")
            self.assertEqual(snp["alt_fragmentos"], "3")
            self.assertEqual(
                snp["confirmacao_bam"],
                "SUPORTE_COMPATIVEL_COM_GENOTIPO"
            )

            sem_cobertura = por_posicao[pos_sem_cobertura]
            self.assertEqual(
                sem_cobertura["status"],
                "SEM_COBERTURA_INICIAL"
            )
            self.assertEqual(
                sem_cobertura["confirmacao_bam"],
                "NAO_AVALIAVEL_SEM_COBERTURA"
            )

            insercao = por_posicao[pos_insercao]
            self.assertEqual(
                insercao["tipo"],
                "INSERCAO_SIMPLES"
            )
            self.assertEqual(
                insercao["status"],
                "TIPO_NAO_SUPORTADO"
            )
            self.assertEqual(
                insercao["confirmacao_bam"],
                "NAO_AVALIAVEL_TIPO"
            )


if __name__ == "__main__":
    unittest.main(verbosity=2)
