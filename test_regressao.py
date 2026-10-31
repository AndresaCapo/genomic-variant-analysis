import tempfile
import unittest
from pathlib import Path

import pysam

from exportar_resultados import classificar_tipo
from suporte_alelico import contar_snp, contar_delecao


def criar_bam(leituras):
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
