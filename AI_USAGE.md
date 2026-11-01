# AI_USAGE — Utilização de Inteligência Artificial

## 1. Ferramenta utilizada

**Ferramentas:** ChatGPT (OpenAI) e Claude (Anthropic).

O ChatGPT foi utilizado como assistente técnico durante o desenvolvimento
da solução, apoiando a investigação exploratória dos dados genômicos,
a definição de critérios analíticos, a implementação em Python,
a elaboração de testes e a documentação.

O processo foi iterativo: as sugestões foram executadas, seus resultados
foram examinados e ajustes foram realizados quando necessários.

O Claude foi utilizado posteriormente como segunda ferramenta de IA
para uma revisão crítica do código, da documentação e do relatório
técnico, com foco na identificação de fragilidades metodológicas,
inconsistências e oportunidades de melhoria.


## 2. Estratégia de investigação e desenvolvimento

Antes da implementação da ferramenta, foi realizada uma exploração dos
arquivos fornecidos, com o objetivo de compreender os dados, identificar
limitações e estabelecer critérios para a análise automatizada.

### 2.1. Caracterização dos arquivos

A investigação começou pela inspeção do VCF, BAM e respectivos índices,
incluindo cabeçalhos, cromossomos, montagem de referência e integridade.

Foram utilizados comandos de samtools, bcftools e verificações MD5.

Foi identificada uma divergência entre o MD5 calculado para o VCF e
o checksum informado no arquivo fornecido, embora o teste de integridade
da compressão tenha sido satisfatório.

Essa etapa orientou a definição dos requisitos de entrada e o registro
das ressalvas de integridade na documentação.

### 2.2. Identificação das variantes nos genes de interesse

Após a identificação da montagem GRCh38 como referência de trabalho,
foram utilizadas coordenadas gênicas provenientes das anotações
do NCBI Gene para FANCA, HEXA e TGFBR2.

As fontes e os intervalos empregados estão registrados no README.md.

As consultas ao VCF identificaram 52 registros nas regiões selecionadas:
36 em FANCA, 7 em HEXA e 9 em TGFBR2.

A investigação também considerou os tipos de variantes, os genótipos
e os campos de profundidade e filtragem.

Essa etapa definiu as regiões de interesse e os campos do CSV final.

### 2.3. Cobertura e suporte alélico no BAM

Em seguida, foi investigada a disponibilidade de alinhamentos
para avaliar as variantes identificadas no VCF.

Apenas quatro das 52 variantes apresentaram cobertura na posição
genômica inicial no BAM filtrado disponibilizado.

Essa observação motivou a distinção entre ausência de cobertura e
ausência de variante, evitando interpretações não sustentadas pelos dados.

Também foram investigadas diferenças entre contagens presentes no VCF
e contagens calculadas diretamente nos alinhamentos do BAM.

Como resultado, a ferramenta passou a registrar explicitamente
o status SEM_COBERTURA_INICIAL.

### 2.4. Qualidade dos alinhamentos e análise de deleções

Foram investigados os efeitos da qualidade de mapeamento (MAPQ),
qualidade das bases (BaseQ), duplicatas e sobreposição de leituras
pareadas nas contagens de suporte aos alelos REF e ALT.

A deleção FANCA chr16:89765073 CAG>C exigiu uma investigação
específica das operações CIGAR, bases de ancoragem e alinhamentos
que não permitiam classificação segura.

Também foi realizada uma investigação exploratória da relação
entre essa deleção e um SNP próximo, a partir de fragmentos
informativos para ambas as posições.

Essas análises orientaram a implementação de funções específicas
para SNPs e deleções simples, com contagens por alinhamento e
por QNAME e identificação de evidências discordantes.

A análise exploratória de fase não integra a saída automatizada.

### 2.5. Implementação da ferramenta

Com os critérios definidos, a investigação passou à elaboração
dos scripts Python, organização dos módulos, desenvolvimento
da interface de linha de comando e geração de resultados em CSV.

O ambiente Docker foi preparado para permitir execução isolada,
com volumes para leitura dos dados e gravação dos resultados.

A ferramenta também recebeu identificação explícita de versão,
consultável pelo argumento --version.

## 3. Principais prompts e perguntas utilizados

Os exemplos a seguir são resumos fiéis dos temas abordados nas
interações, e não transcrições literais de mensagens.

| Etapa | Pergunta ou solicitação resumida | Objetivo |
|---|---|---|
| Inspeção | Como verificar cabeçalhos, índices, montagem genômica e integridade do VCF e BAM? | Caracterizar os arquivos de entrada. |
| Anotação gênica | Como localizar os intervalos de FANCA, HEXA e TGFBR2 na referência GRCh38? | Definir as regiões para filtragem. |
| Variantes | Como identificar e exportar as variantes desses genes usando Python e pysam? | Desenvolver a seleção das variantes. |
| Cobertura | Como verificar quais variantes possuem cobertura no BAM e interpretar posições sem leituras? | Evitar conclusões incorretas. |
| Suporte alélico | Como calcular suporte REF/ALT considerando MAPQ, BaseQ e leituras pareadas? | Definir os critérios de contagem. |
| Deleções | Como identificar deleções pelas operações CIGAR e tratar alinhamentos não classificados? | Implementar suporte para deleções simples. |
| Testes e Docker | Como criar testes sintéticos e executar a ferramenta em ambiente isolado com volumes? | Verificar funcionamento e reprodutibilidade. |
| Versionamento | Como permitir a consulta da versão da aplicação pela linha de comando? | Atender ao requisito de identificação da versão. |

## 4. Validação das respostas e ajustes realizados

As sugestões produzidas pela IA não foram incorporadas automaticamente.
O desenvolvimento incluiu verificações independentes e testes práticos.

### 4.1. Verificação com ferramentas de bioinformática

A seleção de variantes foi conferida com consultas utilizando bcftools.

As evidências nos alinhamentos foram investigadas com samtools,
pysam e inspeção das operações CIGAR.

As diferenças observadas entre profundidades e contagens motivaram
investigações adicionais sobre filtros e leituras sobrepostas.

Como resultado, a ferramenta passou a apresentar separadamente
as contagens por alinhamento e por QNAME.

### 4.2. Ajustes motivados pela investigação

Algumas observações exigiram revisão das abordagens inicialmente
consideradas:

- Cobertura limitada: foram preservados os registros sem cobertura
  inicial, sem classificá-los como variantes ausentes.
- Deleção em FANCA: foram incorporadas categorias para alinhamentos
  não classificados e fragmentos discordantes.
- Leituras pareadas: as contagens por alinhamento e fragmento foram
  diferenciadas para evitar interpretações inadequadas.
- Qualidade: foram definidos filtros de MAPQ e BaseQ e exclusões
  de categorias específicas de alinhamentos.
- FILTER do VCF: variantes com filtros diferentes de PASS foram
  preservadas e identificadas na saída.

Esses ajustes decorreram da análise dos resultados e não apenas
da aceitação das primeiras sugestões da IA.

### 4.3. Testes automatizados

A suíte foi ampliada para 21 testes automatizados, utilizando
dados sintéticos:

- Dez testes iniciais de classificação de variantes, suporte
  alélico, qualidade e interpretação de SNPs e deleções.
- Dez testes para os critérios operacionais de compatibilidade
  do suporte no BAM com o genótipo informado no VCF.
- Um teste de integração que cria VCF e BAM sintéticos, executa
  a exportação completa e verifica o CSV resultante.

Os 21 testes foram executados com sucesso no Ubuntu e no Docker.

A execução completa com os arquivos fornecidos também confirmou
a preservação das 20 colunas anteriores do CSV, acrescentando
somente a coluna `confirmacao_bam`.

### 4.4. Reprodutibilidade

A ferramenta foi executada localmente e em ambiente Docker.

Os resultados CSV gerados foram comparados byte a byte e
permaneceram idênticos para o conjunto de entrada analisado.

A versão 1.4.0 foi identificada e testada por meio de --version,
incluindo sua execução no Docker.

### 4.5. Revisão cruzada com Claude

Após a implementação inicial, o projeto foi submetido a uma
segunda avaliação assistida por IA, utilizando Claude (Anthropic).

As solicitações incluíram:

- Revisar criticamente a seleção de variantes e a análise do BAM.
- Identificar possíveis erros na contagem de SNPs e deleções.
- Avaliar a coerência entre profundidade e suporte alélico.
- Examinar os testes, o Docker, o README e o relatório técnico.
- Verificar o atendimento aos requisitos do desafio.

A revisão destacou a necessidade de explicar as diferenças entre
DP BAM e contagens REF/ALT, definir um critério operacional de
compatibilidade do suporte com o genótipo e ampliar os testes
automatizados.

Como resultado, foi implementada a coluna `confirmacao_bam`,
foi acrescentada a explicação metodológica sobre profundidade
e foi desenvolvido um teste de integração com dados sintéticos.

As alterações foram verificadas posteriormente por testes no
Ubuntu e no Docker e por comparação com os resultados anteriores.

A avaliação realizada com Claude constituiu uma revisão cruzada
assistida por outra IA, não uma auditoria humana independente.
As recomendações foram analisadas criticamente antes de sua
incorporação ao projeto.

## 5. Limitações e responsabilidade

A inteligência artificial foi utilizada como apoio técnico,
mas suas respostas não constituem validação independente.

A ferramenta não recalcula genótipos, não classifica patogenicidade
e não realiza interpretação clínica.

A análise automática de suporte alélico contempla SNPs bialélicos
e deleções simples bialélicas compatíveis com o algoritmo.

As contagens por QNAME não representam necessariamente moléculas
independentes, e a cobertura limitada do BAM restringe as conclusões.

Os testes realizados demonstram o funcionamento nos cenários
avaliados, mas não substituem validação laboratorial independente.

Os arquivos genéticos originais não foram incorporados à imagem
Docker nem ao repositório Git.

## 6. Considerações finais

A IA foi utilizada ao longo de um processo de investigação,
implementação, teste e documentação.

As perguntas iniciais buscaram compreender os arquivos e as
evidências disponíveis antes de iniciar a programação.

Posteriormente, as observações obtidas orientaram decisões
de implementação e verificações de reprodutibilidade.

Essa abordagem permitiu utilizar a IA como apoio ao raciocínio
técnico, mantendo a necessidade de verificar criticamente
as sugestões e reconhecer as limitações dos resultados.
