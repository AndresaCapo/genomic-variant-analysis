# Desafio de Bioinformática — Análise de Variantes Genéticas

## 1. Objetivo

Ferramenta de linha de comando desenvolvida em Python 3 para identificar
variantes nos genes FANCA, HEXA e TGFBR2 e avaliar o suporte alélico
nos alinhamentos disponíveis em um arquivo BAM.

A ferramenta gera um arquivo CSV contendo informações do VCF,
cobertura no BAM e contagens de suporte aos alelos REF e ALT.

## 2. Tecnologias

- Python 3.11
- pysam 0.23.0
- Docker
- unittest

## 3. Arquivos do projeto

- analise_variantes.py: consulta às regiões genômicas e cobertura.
- suporte_alelico.py: contagem de SNPs e deleções simples.
- exportar_resultados.py: geração do arquivo CSV.
- test_regressao.py: testes automatizados.
- Dockerfile: definição do ambiente de execução.
- .dockerignore: exclusão dos dados genéticos da imagem.

## 4. Referência genômica

Montagem: GRCh38/hg38, conforme cabeçalho do VCF.

Regiões utilizadas:

- TGFBR2: chr3:30606356-30694142
- HEXA: chr15:72340924-72376014
- FANCA: chr16:89737549-89816647

As coordenadas estão fixadas no código e devem ser revisadas
caso seja utilizada outra versão da anotação gênica.

## 5. Entradas

O programa recebe:

- VCF compactado em BGZF, com índice CSI ou TBI.
- BAM ordenado por coordenadas, com índice BAI.

Os arquivos devem utilizar a mesma montagem genômica e
nomes de cromossomos compatíveis.

## 6. Construção da imagem Docker

Na pasta do projeto:

    sudo docker build -t desafio-dasa-bioinfo:1.3 .

O arquivo .dockerignore evita incluir os dados genéticos
no contexto de construção.

## 7. Execução

Criar a pasta de saída:

    mkdir -p saida_docker

Executar:

    sudo docker run --rm \
      --network none \
      --user "$(id -u):$(id -g)" \
      --mount type=bind,source="$(pwd)",target=/dados,readonly \
      --mount type=bind,source="$(pwd)/saida_docker",target=/saida \
      desafio-dasa-bioinfo:1.3 \
      --vcf /dados/desafio.vcf.gz \
      --bam /dados/desafio.bam \
      --saida /saida/resultados_variantes.csv

Os arquivos de entrada são montados somente para leitura.
A execução ocorre sem acesso à rede.

## 8. Metodologia

As variantes são selecionadas pelas posições iniciais
situadas nas regiões genômicas definidas para os três genes.

O programa preserva informações originais do VCF,
incluindo genótipo, FILTER, DP e AD.

A cobertura inicial é calculada separadamente do suporte alélico.

A avaliação de suporte utiliza:

- MAPQ mínimo: 20.
- BaseQ mínimo: 20.
- Exclusão de alinhamentos duplicados, secundários,
  suplementares, não mapeados ou com falha de QC.
- Identificação de bases REF e ALT para SNPs bialélicos.
- Inspeção de bases e operações CIGAR para deleções
  simples bialélicas com base de ancoragem inicial.
- Contagens separadas por alinhamento e por QNAME.
- Identificação de evidências discordantes por fragmento.

As contagens são baseadas nas qualidades originais das
leituras e não reproduzem necessariamente o tratamento
de sobreposição utilizado pelo samtools mpileup.

## 9. Interpretação dos resultados

SUPORTE_AVALIADO:
O programa conseguiu aplicar a análise alélica implementada.

SEM_COBERTURA_INICIAL:
Não foi identificada cobertura na posição inicial da variante.

TIPO_NAO_SUPORTADO:
A variante não possui um formato contemplado pelo
algoritmo atual de suporte alélico.

A ausência de cobertura não demonstra ausência da variante.

FILTER=PASS representa o resultado da filtragem original
do VCF, não uma confirmação independente.

## 10. Resultados do conjunto fornecido

Foram identificados 52 registros:

- FANCA: 36
- HEXA: 7
- TGFBR2: 9

Classificação estrutural:

- SNPs: 42
- Deleções simples: 8
- Inserções simples: 2

Cobertura na posição inicial:

- Com cobertura: 4
- Sem cobertura: 48

Das quatro variantes com cobertura, três são SNPs e
uma é uma deleção simples.

## 11. Testes

Executar os testes dentro do Docker:

    sudo docker run --rm \
      --network none \
      --entrypoint python3 \
      desafio-dasa-bioinfo:1.3 \
      /app/test_regressao.py

Foram aprovados dez testes automatizados, incluindo verificacoes de qualidade, fragmentos discordantes e delecoes com soft clipping.

A execução local e a execução Docker produziram
arquivos CSV idênticos no conjunto utilizado.

## 12. Limitações

- A avaliação de cobertura é baseada na posição inicial.
- A ausência de cobertura nessa posição não exclui
  cobertura parcial em outras partes de variantes complexas.
- O suporte alélico não está implementado para todos
  os tipos possíveis de variantes.
- A contagem por QNAME não representa necessariamente
  o número de moléculas independentes.
- Fragmentos discordantes exigem interpretação cautelosa.
- As chamadas do VCF não são recalculadas pelo programa.
- A ferramenta não realiza classificação de patogenicidade.
- O programa não executa interpretação clínica.
- O faseamento local investigado exploratoriamente
  não integra a exportação automatizada.
- A análise não substitui validação laboratorial independente.

## 13. Integridade e privacidade

Os arquivos genéticos não são incorporados à imagem Docker.

O processamento pode ser realizado inteiramente em ambiente
local, sem acesso à internet.

Foi observada uma divergência na verificação MD5 do VCF
em relação ao arquivo md5sum.txt fornecido. A integridade
de compressão foi verificada separadamente. A divergência
de checksum deve ser esclarecida antes de considerar
a integridade do conjunto original completamente confirmada.

## 14. Uso de inteligência artificial

O desenvolvimento contou com assistência de inteligência
artificial para elaboração, revisão e testes do código.

As atividades e limitações dessa assistência serão
documentadas separadamente no arquivo AI_USAGE.md.
