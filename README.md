# # Desafio de Bioinformática — Análise de Variantes Genéticas

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

### Origem das coordenadas genômicas

Após a identificação da montagem GRCh38 nos arquivos
fornecidos, foram consultadas as anotações do NCBI Gene
(RefSeq RS_2025_08, GRCh38.p14) para definir os intervalos
genômicos dos três genes de interesse.

Referências:

- TGFBR2: NCBI Gene ID 7048
  https://www.ncbi.nlm.nih.gov/gene/7048
  chr3:30606356-30694142

- HEXA: NCBI Gene ID 3073
  https://www.ncbi.nlm.nih.gov/gene/3073
  chr15:72340924-72376014

- FANCA: NCBI Gene ID 2175
  https://www.ncbi.nlm.nih.gov/gene/2175
  chr16:89737549-89816647

As coordenadas representam intervalos de anotacao genica
na referencia GRCh38.p14, e nao posicoes obtidas do VCF.

As variantes foram identificadas separadamente no VCF
fornecido pela DASA, utilizando esses intervalos.

## 5. Entradas

O programa recebe:

- VCF compactado em BGZF, com índice CSI ou TBI.
- BAM ordenado por coordenadas, com índice BAI.

Os arquivos devem utilizar a mesma montagem genômica e
nomes de cromossomos compatíveis.

## 6. Construção da imagem Docker

Na pasta do projeto:

    sudo docker build -t desafio-dasa-bioinfo:1.4.0 .

O arquivo .dockerignore evita incluir os dados genéticos
no contexto de construção.

## 7. Execução

### Identificação da versão

Para consultar a versão da aplicação:

    sudo docker run --rm --network none desafio-dasa-bioinfo:1.4.0 --version

Resultado esperado:

    exportar_resultados.py 1.4.0

Criar a pasta de saída:

    mkdir -p saida_docker

Executar:

    sudo docker run --rm \
      --network none \
      --user "$(id -u):$(id -g)" \
      --mount type=bind,source="$(pwd)",target=/dados,readonly \
      --mount type=bind,source="$(pwd)/saida_docker",target=/saida \
      desafio-dasa-bioinfo:1.4.0 \
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

### Classificacao do suporte computacional no BAM

O CSV inclui a coluna `confirmacao_bam`, que apresenta uma avaliacao
operacional da compatibilidade entre o suporte observado no BAM e
o genotipo registrado no VCF.

Para a classificacao `SUPORTE_COMPATIVEL_COM_GENOTIPO`, sao exigidos:

- FILTER=PASS no VCF.
- Genotipo heterozigoto (0/1, 1/0 ou equivalente faseado):
  pelo menos 3 QNAMEs REF e 3 QNAMEs ALT.
- Genotipo homozigoto alternativo (1/1):
  pelo menos 3 QNAMEs ALT e proporcao ALT de pelo menos 90%
  entre os fragmentos classificados como REF ou ALT.

Outras classificacoes possiveis:

- `SUPORTE_INSUFICIENTE`: evidencia abaixo dos criterios definidos.
- `NAO_AVALIAVEL_SEM_COBERTURA`: ausencia de cobertura inicial no BAM.
- `NAO_AVALIAVEL_TIPO`: tipo de variante sem avaliacao implementada.
- `NAO_AVALIAVEL_FILTRO_VCF`: registro com FILTER diferente de PASS.
- `NAO_AVALIAVEL_GENOTIPO`: genotipo fora do escopo da regra.

Os criterios acima foram definidos operacionalmente para este desafio.
Nao constituem limiares clinicos validados. A contagem por QNAME
nao garante independencia molecular, e a classificacao nao equivale
a confirmacao clinica ou laboratorial independente.

A profundidade inicial no BAM e as contagens REF/ALT utilizam filtros
diferentes; portanto, seus valores nao precisam coincidir.

## 11. Testes

Executar os testes dentro do Docker:

    sudo docker run --rm \
      --network none \
      --entrypoint python3 \
      desafio-dasa-bioinfo:1.4.0 \
      /app/test_regressao.py

Foram aprovados 21 testes automatizados no Ubuntu e no Docker. A suite contempla classificacao de variantes, suporte alelico, qualidade dos alinhamentos, fragmentos discordantes, delecoes, criterios de confirmacao computacional e um teste de integracao com VCF e BAM sinteticos, incluindo a exportacao completa do CSV.

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

As atividades e limitações dessa assistência estão
documentadas separadamente no arquivo AI_USAGE.md.

## 15. Fundamentação metodológica e validação

### Qualidade dos alinhamentos: MAPQ e BaseQ

Foram adotados MAPQ >= 20 e BaseQ >= 20 para a avaliação do
suporte alélico. Na escala Phred, Q20 corresponde a uma
probabilidade estimada de erro de 1%.

MAPQ representa a confiança no posicionamento do alinhamento,
enquanto BaseQ representa a qualidade estimada da chamada
de uma base. Esses limiares são critérios operacionais deste
projeto e não constituem limites clínicos universalmente validados.

### Profundidade de cobertura e suporte alélico

A profundidade inicial no BAM e as contagens REF/ALT são
medidas distintas. A cobertura é calculada na posição inicial
da variante, enquanto o suporte alélico considera filtros
adicionais de qualidade e classificação dos alelos.

Consequentemente, REF + ALT pode ser diferente de DP BAM.

### Identificação de deleções pelo CIGAR

O CIGAR descreve as operações de alinhamento de uma leitura
contra a referência. A operação D indica bases presentes
na referência e ausentes na leitura.

Para deleções simples, o algoritmo verifica a base de
ancoragem e a posição e extensão esperadas da deleção.
Representações alternativas podem não ser reconhecidas.

### Agrupamento por QNAME

Os alinhamentos são agrupados pelo identificador QNAME para
reduzir a contagem duplicada de evidências provenientes de
leituras pareadas sobrepostas. Evidências incompatíveis
podem ser classificadas como discordantes.

QNAME não corresponde necessariamente a moléculas independentes.

### Reprodutibilidade computacional

O repositório foi clonado em um diretório independente e a
imagem Docker foi construída novamente sem utilizar cache.

Os 21 testes automatizados foram aprovados. O processamento
dos arquivos VCF e BAM produziu 52 registros de variantes.

A comparação dos arquivos CSV demonstrou identidade byte a byte
entre o resultado reproduzido e o obtido anteriormente.

### Limites da validação

Os testes realizados demonstram o funcionamento computacional
nos cenários avaliados, mas não estabelecem validade clínica.

O uso diagnóstico exigiria validação adicional com materiais
de referência, conjuntos independentes, avaliação de
sensibilidade, especificidade, precisão e reprodutibilidade,
além de critérios formais de controle de qualidade.
