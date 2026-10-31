# AI_USAGE — Declaração de Uso de Inteligência Artificial

## 1. Ferramenta utilizada

ChatGPT, da OpenAI.

A ferramenta foi utilizada como assistente durante o desenvolvimento
da solução para o desafio de bioinformática.

## 2. Finalidade do uso

A inteligência artificial auxiliou nas seguintes atividades:

- Planejamento das etapas de investigação dos arquivos VCF e BAM.
- Sugestão de comandos para inspeção e controle de qualidade.
- Elaboração de scripts em Python 3.
- Implementação das funções de contagem de suporte alélico.
- Desenvolvimento da exportação dos resultados para CSV.
- Preparação do Dockerfile e do arquivo .dockerignore.
- Elaboração de testes automatizados com dados sintéticos.
- Organização da documentação técnica.
- Discussão das limitações metodológicas da análise.

## 3. Processo de desenvolvimento

O desenvolvimento foi realizado de maneira iterativa.

As instruções e os trechos de código sugeridos pela IA foram
executados no ambiente Ubuntu utilizado no desafio.

Os resultados apresentados pelas ferramentas foram analisados
ao longo do desenvolvimento, permitindo identificar diferenças
entre contagens, limitações de cobertura e situações que
necessitavam de investigação adicional.

Os testes incluíram verificações com ferramentas como
bcftools, samtools, Python, pysam e Docker.

## 4. Verificações realizadas

Foram realizadas as seguintes verificações:

- Inspeção dos cabeçalhos e das estruturas do VCF e do BAM.
- Consulta das variantes nas regiões dos genes selecionados.
- Comparação de contagens alélicas entre VCF e BAM.
- Inspeção de operações CIGAR para uma deleção em FANCA.
- Avaliação exploratória de orientação dos alinhamentos.
- Análise de leituras pareadas e fragmentos discordantes.
- Testes unitários com dados sintéticos.
- Execução dos testes automatizados em Docker.
- Comparação byte a byte dos CSVs produzidos localmente
  e dentro do Docker.

O arquivo permanente test_regressao.py contém cinco testes
automatizados que foram executados com sucesso.

## 5. Resultados obtidos

O programa identificou 52 registros de variantes nas regiões
genômicas utilizadas para FANCA, HEXA e TGFBR2.

Foram identificados:

- 42 SNPs.
- 8 deleções simples.
- 2 inserções simples.

Quatro variantes apresentaram cobertura na posição inicial
e tiveram o suporte alélico avaliado pelo algoritmo implementado.

As demais 48 foram mantidas no resultado como variantes
sem cobertura na posição inicial no BAM disponibilizado.

## 6. Limitações e cuidados

A utilização de inteligência artificial não substitui
verificação experimental ou revisão científica independente.

A interpretação dos resultados depende das definições
adotadas no algoritmo e da cobertura dos alinhamentos.

O código implementado não realiza classificação clínica
de patogenicidade nem confirmação laboratorial de variantes.

As contagens por QNAME não representam necessariamente
moléculas independentes.

A avaliação automática de suporte alélico está limitada
a SNPs bialélicos e deleções simples bialélicas no formato
suportado pelo programa.

O faseamento local explorado durante o desenvolvimento
não integra a saída automatizada.

As coordenadas gênicas e os critérios de análise devem
ser considerados explicitamente ao interpretar os resultados.

Foi observada divergência entre o MD5 do VCF e o checksum
informado no arquivo md5sum.txt fornecido.

Os testes realizados demonstram funcionamento nos cenários
avaliados, mas não garantem ausência de erros em todos
os casos possíveis.

## 7. Privacidade

O desenvolvimento foi realizado utilizando ferramentas
locais de análise.

Os arquivos genéticos não foram incorporados à imagem Docker.

A execução do programa no contêiner foi realizada sem
acesso à rede e com os dados de entrada montados
somente para leitura.

Os testes automatizados utilizam dados sintéticos.

## 8. Responsabilidade e transparência

Os códigos e documentos receberam assistência de IA.

A execução dos comandos, a conferência dos resultados
e a preparação do ambiente foram realizadas durante
o desenvolvimento do desafio.

A presente declaração tem como finalidade registrar
essa assistência de forma transparente.

A reprodutibilidade computacional observada não deve
ser confundida com validação clínica ou laboratorial.
